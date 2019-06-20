#!/bin/bash

EID=''
SALT=''

function help() {
  echo "usage: ehash [options]"
  echo -e "  -i Employee ID (decimal integer)"
  echo -e "  -s Salt (HEX number)"
  echo -e "  -h Show this message"
  exit 0
}

while getopts 'i:s:h' flag; do
  case "${flag}" in
    i) EID="${OPTARG}" ;;
    s) SALT="${OPTARG}" ;;
    h) help ;;
    *) help ;;
  esac
done

if [[ -z "${EID}" ]] || [[ -z ${SALT} ]]; then
  echo "Arguments missing..."
  exit 1
fi

# Transform arguments to binary numbers
SALT=$(echo "ibase=16;obase=2;${SALT}" | bc | tr -d '\\\n')
EID=$(echo "ibase=10;obase=2;${EID}" | bc | tr -d '\\\n')

# Concat two bit strings but first pad them to make each 256 bits long
PSALT=$(printf "%0256s" "${SALT}" | tr ' ' '0')
PEID=$(printf "%0256s" "${EID}" | tr ' ' '0')
PAYLOAD="${PEID}${PSALT}"

# Print only hash (and not '-' for stdout)
echo $PAYLOAD | shasum -a 256 --01 | cut -d ' ' -f 1

