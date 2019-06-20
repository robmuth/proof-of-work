#!/bin/bash

EID=''
# Get script directory
# https://stackoverflow.com/a/246128/2295964
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

function help() {
  echo "usage: stamp [options]"
  echo -e "  -i Employee ID (decimal integer)"
  echo -e "  -h Show this message"
  exit 0
}

while getopts 'i:h' flag; do
  case "${flag}" in
    i) EID="${OPTARG}" ;;
    h) help ;;
    *) help ;;
  esac
done

if [[ -z "${EID}" ]]; then
  echo "Arguments employee ID..."
  exit 1
fi

SALT=$($DIR/salt.sh)
HASH=$($DIR/ehash.sh -i $EID -s $SALT)
TERMSTAMP=$($DIR/terminal.py "${HASH}")
PSALT=$(echo "ibase=16;obase=2;${SALT}" | bc | tr -d '\\\n' | sed 's/./& /g')

echo "${TERMSTAMP} ${EID} ${PSALT%% }"