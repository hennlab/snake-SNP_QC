REPORT=$1
DATASET=$2
OUTPUT=$3

BESTZ=`grep -oP '(?<=Global Concordance: ).*' zCall/*concordance.stats* | sort -sn -t ":" -k 2,2 | tail -n 1 | awk -F ":" '{print$1}' | sed -n 's/.*\(stats[0-9]\).*/\1/p' | sed -e s/[^0-9]//g`
echo $BESTZ
python scripts/zCall.py -R ${REPORT} -T zCall/${DATASET}.thresholds.${BESTZ}.txt -O ${OUTPUT}
