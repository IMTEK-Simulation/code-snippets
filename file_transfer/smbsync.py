#!/usr/bin/env python3
""" One-way sync to an SMB/CIFS server (needs pysmb)"""
# 9 Aug 2018, adapted for Uni Freiburg's RZ Storage "pastewka" share
# No NetBIOS needed, direct SMB via TCP 445
# Compares files only based on file size.
# TODO:
#  * Combine 'standard' and modified version DONE, Jun 2019
#  * Add some 'protected' password method, avoid command line argument DONE, Jun 2019

import logging, socket, os, sys, math, configparser
from datetime import datetime

import logging
logger = logging.getLogger(__name__)
logfmt = "[%(levelname)s - %(filename)s:%(lineno)s - %(funcName)s() ] %(message)s (%(asctime)s)"
logging.basicConfig( format = logfmt )

try:
    from smb.SMBConnection import SMBConnection
    from smb.base import NotReadyError, NotConnectedError
    from nmb.NetBIOS import NetBIOS
    PYSMB = True
except ImportError:
    PYSMB = False

def props_SharedFile(sf):
    return {
      'ctime': datetime.utcfromtimestamp(sf.create_time),
      'mtime': datetime.utcfromtimestamp(sf.last_write_time),
      'size':  math.ceil(sf.file_size/1024.0),
      'file':  not sf.isDirectory,
      'name':  sf.filename,
    }

def props_LocalFile(path):
    return {
      'mtime': datetime.utcfromtimestamp(os.path.getmtime(path)),
      'ctime': datetime.utcfromtimestamp(os.path.getctime(path)),
      'size':  math.ceil(os.path.getsize(path)/1024.0),
      'file':  os.path.isfile(path),
      'name':  os.path.basename(path),
    }

def info_for_SharedFile(sf):
    props = props_SharedFile(sf)
    props['kind'] = 'file' if props['file'] else ' dir'
    tpl = "Created: {ctime}, Last change: {mtime}, Size: {size} KiB, Kind: {kind}, Name: {name}"
    return tpl.format(**props)

def smbsync(
    from_path, to_path,
    credentials = os.path.expanduser('~/.smbcredential.rz_storage'),
    server_ip   = None,
    share       = "pastewka",
    domain      = "PUBLIC",
    username    = "pastewka",
    password    = None,
    server_name = "ufr2.isi1.public.ads.uni-freiburg.de",
    server_port = 445 ):

    # read credentials from file if specified and existant:
    if os.path.isfile(credentials):
        # create dumme section, otherwise ConfigParser throws exception:
        config_str = '[root]\n' + open(credentials,'r').read()
        logger.info( "Found credentials at '{:s}', parsing...".format(credentials) )
        config = configparser.ConfigParser()
        config.read_string(config_str)
        if 'username' in config['root']:
            logger.debug( "Found username '{:s}' in credentials file.".format(
              config['root']['username'] ) )
            if username is not None:
              logger.debug( "Overriding previously specified username '{:s}'.".format(
                username ) )
            username = config['root']['username']
        if 'password' in config['root']:
            logger.debug( "Found password '{:s}' in credentials file.".format(
              config['root']['password'] ) )
            if password is not None:
              logger.debug( "Overriding previously specified password '{:s}'.".format(
                password ) )
            password = config['root']['password']
    try:
        #nb = NetBIOS(broadcast=True)
        #if not server_name:
        #    server_name = nb.queryIPForName(server_ip)[0]
        if not server_ip:
            server_ip = socket.gethostbyname(server_name)
        #    server_ip = nb.queryName(server_name)

        # conn = SMBConnection(username, password, socket.gethostname(), server_name, use_ntlm_v2=True, domain=domain)
        host_name = socket.gethostname()
        conn = SMBConnection(username, password, host_name,
                server_name, domain=domain, 
                use_ntlm_v2=True, is_direct_tcp=True)

        logger.info( ( "Connecting from '{host:s}' to "
            "'smb://{user:s}@{ip:s}({server:s}):{port:d}', "
            " DOMAIN '{domain:s}'").format(user=username,
                ip=server_ip, server=server_name, 
                port=server_port, host=host_name, 
                domain=domain) )

        # for testing, see types of arguments
        logger.debug( ( "Types HOST '{host:s}', USER '{user:s}', IP '{ip:s}', "
           "SERVER '{server:s}', PORT '{port:s}', DOMAIN '{domain:s}', "
            "PASSWORD '{password:s}'").format(
                user=type(username).__name__, 
                ip=type(server_ip).__name__, 
                server=type(server_name).__name__, 
                port=type(server_port).__name__, 
                host=type(host_name).__name__, 
                domain=type(domain).__name__,
                password=type(password).__name__ ) )

        # conn.connect(server_ip, 139, timeout=30)
        conn.connect(server_ip, port=server_port, timeout=30)

    except (NotReadyError, NotConnectedError) as e:
        logger.error('Could not connect: ' + e.__class__.__name__)
        return 2
    except Exception as e:
        logger.error('Could not connect:')
        logger.error(e.__class__.__name__ + " " + str(e))
        return 2
    logger.info("Connection established.")

    leading_path = os.path.normpath(to_path).split(os.sep)
    trailing_path = []
    while True:
        try: # first check whether directory exists ...
            filelist = conn.listPath( share, os.path.join(*leading_path) )
            logger.info("Remote directory '{:s}' exists.".format(os.path.join(*leading_path)))
            if len(trailing_path) == 0:
                break # all nested direcories exist
            leading_path.append( trailing_path.pop(0) )
        except: # directory does not exist:
            try: # try to create
                logger.info("Remote directory '{:s}' does not exist, trying to create.".format(os.path.join(*leading_path)))
                conn.createDirectory(share, os.path.join(*leading_path))
            except: # directory could not be created, go up in hierarchy
                logger.info("Failed to create remote directory '{:s}'.".format(os.path.join(*leading_path)))
                trailing_path.insert(0, leading_path.pop())
                if len(leading_path) == 0: # arrived at root of share and did not succeed in creating any directory
                    logger.error("Desired directory '{:s}' could not be created.".format(to_path))
                    return 2

    filelist = [f for f in filelist if f.filename not in ['.', '..']]
    filename_index = dict()
    logger.info("Files currently stored on the server:")
    for sf in filelist:
        logger.info(info_for_SharedFile(sf))
        filename_index[sf.filename] = sf
    # Loop over local files
    logger.info("Files in local folder to sync:")
    for fn in os.listdir(from_path):
        if not os.path.isfile(os.path.join(from_path, fn)):
            logger.info('Skipping local folder: {}'.format(fn))
            continue
        if fn not in filename_index:
            logger.info("{} not yet on the server. uploading...".format(fn))
            conn.storeFile(share, os.path.join(to_path, fn), open(os.path.join(from_path, fn), 'rb'))
        elif props_SharedFile(filename_index[fn])['size'] != props_LocalFile(os.path.join(from_path, fn))['size']:
            logger.info("{} outdated or incomplete on the server. updating...".format(fn))
            conn.storeFile(share, os.path.join(to_path, fn), open(os.path.join(from_path, fn), 'rb'))
        else:
            logger.info("{} already on the server.".format(fn))

    return 0

def main():
    if not PYSMB:
        print('Could not find the dependency pysmb. Please install it first using `pip install pysmb`.')
        sys.exit(127)

    import argparse
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--credentials', '-c', help='Config file to read credentials (username and password) from',
        default=os.path.expanduser('~/.smbcredentials'))
    parser.add_argument('--server_ip', '-i', help='SMB/CIFS Server IP (optional)')
    parser.add_argument('--share', '-s', help='Name of the share"',
        default="pastewka")
    parser.add_argument('--domain', '-d',
        help='Often "WORKGROUP" or "PUBLIC"', default="WORKGROUP")
    parser.add_argument('--username', '-u',
        default="smbuser")
    parser.add_argument('--password', '-p')
    parser.add_argument('--server_name', '-n',
        help='SMB/CIFS Server Name',
        default="localhost")
    parser.add_argument('--server_port', help='Port number, default 445',
        type=int, default=445)
    parser.add_argument('--verbose', '-v', action='store_true',
        help='Make this tool more verbose')
    parser.add_argument('--debug', action='store_true',
        help='Make this tool print debug info')
    parser.add_argument('from_path', help='Local path to sync to the share')
    parser.add_argument('to_path', help='Path in the share to sync to')
    args = parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    elif args.verbose:
        loglevel = logging.INFO
    else:
        loglevel = logging.WARNING
    logger.setLevel(loglevel)

    sys.exit( smbsync(
        from_path   = args.from_path,
        to_path     = args.to_path,
        credentials = args.credentials,
        server_ip   = args.server_ip,
        share       = args.share,
        domain      = args.domain,
        username    = args.username,
        password    = args.password,
        server_name = args.server_name,
        server_port = args.server_port ) )

if __name__ == "__main__":
    main()
