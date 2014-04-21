awk 'BEGIN{FS=","}; {for(i=1;i<=NF;i++) { max[i]=($i > max[i] ? $i : max[i]); min[i]=($i < min[i] ? $i : min[i]);}} END {for(i=1;i<=NF;i++) { printf("%f,", min[i]); } print ""; for(i=1;i<=NF;i++) { printf("%f,", max[i]);} print "";}'
#tail ../uk3day/train.csv | cut -d , -f 2 | sort -nr
#tail ../uk3day/train.csv | sh scalepass.sh