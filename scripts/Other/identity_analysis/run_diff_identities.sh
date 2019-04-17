

## a script to run the GROUP task of SLING with varying identity levels so I could compared their results

array=( 35 45 55 65 75 85 95 )

for i in "${array[@]}"
do
	bsub -J ${i}_ident -G team216 -o ${i}.o -e ${i}.e -R"select[mem>5000] rusage[mem=5000]" -M5000 -n6 -R"span[hosts=1]"  sling group newFilter ${i} toxins -u -c 6 -it -mi ${i}
done