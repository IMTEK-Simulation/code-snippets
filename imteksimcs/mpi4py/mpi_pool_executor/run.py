"""Run callable MPI-parallel by fully qualified name."""
import argparse
import logging

from six.moves import builtins
from .mpi_pool_executor import call


class ArgumentDefaultsAndRawDescriptionHelpFormatter(
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.RawDescriptionHelpFormatter):
    """Allows for both preformatted help and automatic defaults display."""
    pass


def resolve_fqn(fqn):
    """Resolve fully qulified name to callable object."""
    toks = fqn.rsplit('.', 1)
    if len(toks) == 2:
        modname, funcname = toks
        mod = __import__(modname, globals(), locals(), [str(funcname)], 0)
        func = getattr(mod, funcname)
    else:
        # Handle built in functions.
        func = getattr(builtins, toks[0])
    return func


def main():
    """Run MPI-parallel callable by fully qualified name from command line."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=ArgumentDefaultsAndRawDescriptionHelpFormatter)

    # global logging options
    parser.add_argument('--debug', '-d', default=False, required=False,
                        action='store_true', dest="debug", help='debug')
    parser.add_argument('--verbose', '-v',  default=False, required=False,
                        action='store_true', dest="verbose", help='verbose')

    parser.add_argument('--log', '-l', required=False, nargs='?', dest="log",
                        default=None, const='out.log', metavar='LOG',
                        help='Write out.log, optionally specify log name')

    # positional options
    parser.add_argument('callable', metavar='CALLABLE',
                       help="Callable (i.e. function by fully qualified name.")

    parser.add_argument('args', nargs=argparse.REMAINDER)

    args = parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    elif args.verbose:
        loglevel = logging.INFO
    else:
        loglevel = logging.WARNING

    logging.basicConfig(level=loglevel)

    # explicitly modify the root logger (necessary?)
    logger = logging.getLogger()
    logger.setLevel(loglevel)

    # log everything to file if desired
    if args.log:
        fh = logging.FileHandler(args.log)
        fh.setLevel(loglevel)
        logger.addHandler(fh)

    logger.info("This is '{}' : '{}'.".format(__file__,__name__))

    for arg, val in sorted(vars(args).items()):
        logging.info("Argument {arg:}: '{val:}' ({type:})".format(arg=arg,val=val,type=type(val)))

    func = resolve_fqn(args.callable)
    call(func, *args.args)