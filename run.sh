A=0
while [ $A -le $1 ]
do
	./modelEval&
	A=`expr $A + 1`
done
