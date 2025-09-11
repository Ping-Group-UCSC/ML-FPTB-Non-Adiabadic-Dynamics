for i in 1 2 3
do
  cd $i 
  /g/g92/calegari/deepmd-kit_from_sh/bin/dp freeze -o graph${i}.pb
  cd ..
done  
