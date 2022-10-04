# Copy all datasets in a dependency graph to another base URI.
import dtoolcore
from dtool_lookup_api.core.LookupClient import ConfigurationBasedLookupClient

cl = ConfigurationBasedLookupClient()

await cl.connect()

# Query dependency graph of dataset identified by 'uuid' from dtool lookup server
uuid = '6d7633f-c629-4d30-ad5b-88a64cf654de'
ret = await cl.graph(uuid)

# extract list of URIs

uri_list = [e["uri"] for e in ret]

# copy all datasets in graph to destination base URI
for uri in uri_list:
    try:
        dtoolcore.copy(uri, 's3://frct-livmats')
    except Exception as e1:
        print(e1)
        try:
            dtoolcore.copy_resume(uri, 's3://frct-livmats')
        except Exception as e2:
            print(e2)
