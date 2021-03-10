for FREQ in $(ls *.nVR_freq)
do
  echo $FREQ
  ./translate-nVR_freq.py $FREQ
done
