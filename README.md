# pneumoSerotyper

NAME

		pneumoSerotyper.py - Tool for identifying pneumococcal serotypes with Nucmer.

SYNOPSIS

		pneumoSerotyper.py -i seq1,seq2,...,seqN -o output

DESCRIPTION

		Run Nucmer on the provided whole genome sequences and process the matches to determine 
		how they similar they are to the pneumococcal capsule polysaccharise synthesis (cps) locus
		sequence which encodes the outer cell polysaccharides that determines the serotypes.

		By default the input file format for the sequences is FASTA.

		The user is required to supply the input files and optionally specify
		other Options. The program can be run as follows;

		#Simplest way to run it is to provide the input sequences in fasta format (default).
		pneumoSerotyper.py -i *.fasta -o output

		#You can optionally specify other options minimum match length to be considered
		pneumoSerotyper.py -i *.fasta -o output -c 500

		A summary of the output file names is given at the end of program's execution.

OPTIONS

		-h
		Help. Print this help information and exit.

		-i filenames
		Specify the input fasta sequence files for the program. You can use wildcard (*)
		to speficy multiple files.

		-o filename
		Specify the prefix for the output files. By default the output files will have
		a prefix "Serotypes.Summary" if not specified by the user.

		-c value
		Specify the minimum length of the matches to be considered (default:500). Note that
		smaller values will result in more spurious matches.

		-r
		Remove Nucmer output files (default is to keep the files)

AUTHOR

		Chrispin Chaguza, Chrispin.Chaguza@liverpool.ac.uk. June 2015

FILES

		pneumoSerotyper.py

DEPENDENCIES

		biopython (http://biopython.org/wiki/Main_Page) 
		http://mummer.sourceforge.net/ 
		
		(if you have homebrew installed, these can be install using commands; 'brew install mummer'
		and 'brew install biopython')
