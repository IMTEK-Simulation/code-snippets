# Compare local dataset repository aganst registered entries on lookup server

Sample usage of `dtool_lookup_compare.py` assuming a local folder `DATASETS` containing some datasets:

```console
$ python dtool_lookup_compare.py -q '{"creator_username":"YOUR_USERNAME"}' DATASETS/

Datasets equal on source and target:

[
    "009a5920-1875-46a2-b2f4-42a8f15d7358",
    "85eecd9d-f330-465e-aa77-63341ec8c8ae",
    "05f26707-bc2f-45b6-96d4-80e2728868f9",
    "a2b7c1e0-cadd-4419-a09a-fd25bcecea6e",
    "5b597394-314f-4dcc-a391-969f8bb182ba",
    "9fa1a3f0-c357-4663-9da1-6128951e0a73"
]

Datasets misssing on target:

[
    "5098ac50-8452-453b-b48d-e4c3e403cf36",
    "ed0e7dcb-aacc-411d-8ce2-6d40dea7a366"
]
```

Run `dtool_lookup_compare.py --help` for detailed usage information.

# Sync local to remote dataset repository (with latter indexed by lookup server)
   
To find local datasets missing on lookup server, identify their local directories, and copy to remote destination, one may use the bash snippet

```bash
SOURCE_BASE_URI="/path/to/local/datasets"
TARGET_BASE_URIT="smb://remote-share"
python dtool_lookup_compare.py --missing-only --terse ${SOURCE_BASE_URI} | tr -d '[", ]' | sed '/^$/d' > uuids.txt
dtool ls -v ${SOURCE_BASE_URI} | grep --no-group-separator -B2 -Ff uuids.txt | sed -n 2~3p | awk '{$1=$1};1' > datasets.txt
for d in $(cat datasets.txt); do
    echo $d
    dtool cp $d ${TARGET_BASE_URI}
done
```

Explanation: 

`dtool_lookup_compare.py --missing-only --terse ${SOURCE_BASE_URI}` prints all UUIDs of local datsets not registered yet at the lookup server as a JSON list.
`tr -d '[", ]' | sed '/^$/d'` strips all JSON syntax and yields UUIDs as plain text with one UUID per line.
`dtool ls -v ${SOURCE_BASE_URI}` lists all local datasets by their name together with their full URI in every second and their UUID in every third line.
`grep --no-group-separator -B2 -Ff uuids.txt` greps all matching UUIDs and the previous two lines.
`sed -n 2~3p` yields only every third line beginning from the second line, i.e. URIs.
`awk '{$1=$1};1'` removes leading whitespaces.

If copy fails with `Dataset already exist`error, there might be partially copied datasets (proto datasets in dtool speak) at the remote location.
The lookup server won't register those. Try to resume unfinished copy operations with `dtool cp --resume`, i.e.

    for d in $(cat datasets.txt); do
        echo $d
        dtool cp --resume $d ${TARGET_BASE_URI}
    done

For large amounts of data, run detached and write output to some log file, i.e. by putting the above in a script `dtool_batch_cp.sh` and evoking

    mkdir -p logs
    nohup bash dtool_batch_cp.sh >> "$(pwd)/logs/$(date +%Y%m%d)_dtool_batch_cp.log" 2>&1 &
