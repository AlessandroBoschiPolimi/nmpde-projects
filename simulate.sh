#!/usr/bin/bash

# TODO: Add MESHDIR INPUTDIR OUTPUTDIR

# COMMAND ARGS
mesh_choice="cube"
neumann_func="cube_pull"
forcing_term=0

# PARALLEL STUFF
threads_num=1


help() {

	cat <<- EOF
	Usage: 
	  ./simulate.sh [FLAGS <arg>]
	
	Options:
	  -M <mesh_type>	Choose a mesh type
	  -N <neumann_func>	Choose a neumann function to execute
	  -n <num_threads>	Number of threads to use
	  -F <0|1>		Choose whether to activate the forcing term (rod_bend)
	
	Mesh Types:
	  cube | rod | cup
	
	Neumann Funtions:
	Cube:
	  cube_pull
	Rod:
	Cup:
	EOF

}

parse_arguments() {
	raw_opts="$@"
	parse_state=0
	for opt in $raw_opts; do
		case $parse_state in
			0)
				case $opt in
					# Help
					-h | --help)
						help
						exit 0;;
					# Select Mesh
					-M)
						parse_state=1
						;;
					# Select Neumann Condition
					-N)
						parse_state=2
						;;
					# Select number of threads
					-n)
						parse_state=3
						;;
					# Select forcing term
					-F)
						forcing_term=1
						;;
					*)
						echo "Not a valid option: $opt"
						exit 1;;
				esac
				;;
			1) # Choosing Mesh
				mesh_choice=$opt
				parse_state=0
				;;
			2) # Choosing Neumann Condition
				neumann_func=$opt
				parse_state=0
				;;
			3) # Choosing number of threads
				threads_num=$opt
				echo "Selected number of threads = $threads_num"
				parse_state=0
				;;
		esac
	done
}

run() {
	args="$mesh_choice $neumann_func $forcing_term"

	echo "Running Simulation..."
	if [ $threads_num -eq 1 ]; then
		echo "./PDE-06 $args"
		cd build/ && ./PDE-06 $args && cd ..;
	else
		echo "mpirun -n $threads_num ./PDE-06 $args"
		cd build/ && mpirun -n $threads_num ./PDE-06 $args && cd ..;
	fi
}


parse_arguments $@
run
