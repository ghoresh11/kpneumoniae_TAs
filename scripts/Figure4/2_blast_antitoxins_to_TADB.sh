

for f in *.fasta
do
	echo $f
	blastp -db TADB_Antitoxins.fa -query $f -outfmt 6 -evalue 0.01 -out ${f}_result
done


