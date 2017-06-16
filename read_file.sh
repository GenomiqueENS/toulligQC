config_file=$HOME/Documents/toulligQC-master/config.txt
docker_config_file=$HOME/Documents/toulligQC-master/docker_config.txt
 
index=0
for line in $(cut -f2 -d = "${config_file}"); 
do
	path[index]=$(readlink -f "${line}")
	index=$(expr $index + 1)
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
	
	else
		echo -e "${config_name_var[i]}=${path[i]}/" >> "${docker_config_file}"

		i=$(expr $i+1)
	fi
done


docker run -ti --rm -v ${path[1]}:${path[1]} -v ${path[2]}:${path[2]} -v ${path[3]}:${path[3]} -v ${path[4]}:${path[4]} -v ${path[5]}:${path[5]} genomicpariscentre/toulligqc
