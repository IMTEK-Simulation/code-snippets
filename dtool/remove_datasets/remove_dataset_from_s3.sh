#! /bin/bash

UUID="$1"
ENDPOINT="${2:-https://s3.bwsfs.uni-freiburg.de:8082}"
BUCKET="${3:-frct-simdata}"
AWS_OPTS="${4}" # specify "--profile name-of-profile" to select credentials other than default
DRY_RUN=True

echo "Removing dataset ${UUID} from bucket ${BUCKET} on endpoint ${ENDPOINT}..."

# get prefix from registration key
aws s3 ${AWS_OPTS} --endpoint=${ENDPOINT} cp s3://${BUCKET}/dtool-${UUID} .
PREFIX=`cat dtool-${UUID}`
# rm -f dtool-${UUID}

# remove dataset objects
for fn in `aws s3 ${AWS_OPTS} --endpoint=${ENDPOINT} ls --recursive s3://${BUCKET}/${PREFIX}${UUID}/ | awk '{ print $4 }'`; do
    if [ -n "${DRY_RUN}" ]; then
        echo "Dry run, would remove s3://${BUCKET}/${fn}"
    else
        echo "Remove s3://${BUCKET}/${fn}"
        if ! aws ${AWS_OPTS} s3 --endpoint=${ENDPOINT} rm s3://${BUCKET}/${fn}; then
            echo "Error removing dataset object s3://${BUCKET}/${fn}. Do you have the correct write permissions?"
            exit 1
        fi
    fi
done

# remove registration key
if [ -n "${DRY_RUN}" ]; then
    echo "Dry run, would remove s3://${BUCKET}/dtool-${UUID}"
else
    echo "Remove s3://${BUCKET}/dtool-${UUID}"
    if ! aws ${AWS_OPTS} s3 --endpoint=${ENDPOINT} rm s3://${BUCKET}/dtool-${UUID}; then
        echo "Error removing registration key s3://${BUCKET}/dtool-${UUID}. Do you have the correct write permissions?"
    fi
fi
