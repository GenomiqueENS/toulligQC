config_file=/import/pontos06/analyses/toulligQC/config.txt
docker_config_file=/import/pontos06/analyses/toulligQC/config.txt
 
index=0

for line in $(cut -f2 -d = "${config_file}");

do
  if [ $index -eq 6 ] || [ $index -eq 7 ] || [ $index -eq 8 ]; then
    path[index]=$line
    index=$(expr $index + 1)

  elif [[ $line =~ config ]];then
    path[index]=$line
    index=$(expr $index + 1)

  else
    path[index]=$(readlink -f "${line}")
    index=$(expr $index + 1)
  fi
done

mapfile -t config_name_var < <( cut -f1 -d = ${config_file} )

i=0
if test -f "${docker_config_file}"
       then
               rm "${docker_config_file}"
       fi

for line in ${config_name_var[@]};
do
  if [[ $line =~ config ]]; then
    echo -e [config] >> "${docker_config_file}"
    i=$(expr $i+1)
# echo $line
 elif [[ $line =~ \[extension\] ]]; then
   echo -e [extension] >> "${docker_config_file}"
  i=$(expr $i+1)
  else
    echo -e "${config_name_var[i]}=${path[i]}" >> "${docker_config_file}"
    i=$(expr $i+1)
  fi
done

docker run -ti --rm -v ${path[1]}:${path[1]} -v ${path[2]}:${path[2]} -v ${path[3]}:${path[3]} -v ${path[4]}:${path[4]} -v $(pwd):$(pwd) -v ${path[5]}:${path[5]} -u $(id -u):$(id -g) genomicpariscentre/toulligqc 

