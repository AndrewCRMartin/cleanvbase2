files="allIGHV.fasta allIGKV.fasta allIGLV.fasta allall.fasta humanIGHV.fasta humanIGKV.fasta humanIGLV.fasta humanall.fasta mouseIGHV.fasta mouseIGKV.fasta mouseIGLV.fasta mouseall.fasta"

for file in $files
do
   wget http://www.vbase2.org/download/$file
done
