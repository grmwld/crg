#! /usr/bin/env gawk -f 
#if REF_FILE is defined, filename start_line and end_line are printed in each line. if SCORE is defined, score is printed. If ALI is defined, the alignment is printed also. If ALL is defined, all these are printed. If QLENGTH is defined, the query length is printed after the identity flag
BEGIN {start_pos=-1
q_start_pos=-1
line_start=-1
COMMENT=""
if (ALL){REF_FILE=1
	ALI=1
	SCORE=1
}

}     
{
 
  if ( /^Query=/)   {
	query = $2

	if (line_start==-1)
	{line_start = NR }


  }
  if ( /\(.*letters/)   {

	split($0, LINE_SPLIT, "letters")
	
	while (sub(/[\(\)QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm=\:\ ]/, "", LINE_SPLIT[1])){}

	query_length=LINE_SPLIT[1]



  }

 




if ( /^Sbjct:/ )    
{  
	if (ALI){
		tline=$0
		sub(/Sbjct:/, "")
		while (sub(/[0-9 ]/, "")) {}
		subj_seq = subj_seq $0
		$0=tline 
	}

	while  ( sub(/[A-Z]|[a-z]|\*|\-|:/, " ") )	#    sub(/[:alpha:]|[:upper:]|:/, " ")   )
	{ 	 }
	if ( start_pos==-1 ){	start_pos=$1	}
	end_pos=$2

}
if ( /^Query:/ )    #:[:blank:]*[:digit:]+[:blank:]*[:alpha:]/)
{  
	if (ALI){
		tline=$0
		sub(/Query:/, "")
		while (sub(/[0-9 ]/, "")) {}
		query_seq = query_seq $0
		$0=tline 
	}

	while  ( sub(/[A-Z]|[a-z]|\*|\-|:/, " ") )	#    sub(/[:alpha:]|[:upper:]|:/, " ")   )
	{ 	 }
	if ( q_start_pos==-1 ){	q_start_pos=$1	}
	q_end_pos=$2
}

if (/^>/ || /^\#\#/ || /^ Score/ || /^Reference:/ || /^Matrix:/)
{	
	if (q_start_pos != -1)	#already parsed one hit?
	{ stop_parsing=1}



}
if (stop_parsing==1){
	if (ALI){
		COMMENT=COMMENT " " subj_seq " " query_seq 
		
	}
	if (SCORE)
	{	COMMENT=COMMENT " " evalue " " identity 
	}

	if (QLENGTH)
	{	COMMENT=COMMENT " " query_length 
	}


	if (REF_FILE)
	{	COMMENT=COMMENT " " ARGV[1] " " line_start " " (NR-1) 

	}

	strand="+"
	if (start_pos>end_pos){ temp=start_pos; start_pos=end_pos; end_pos=temp; strand="-"} #exchanging start with stop. using n_hits as a temp var

	print subject " " start_pos " " end_pos " " strand " "  query " " q_start_pos " " q_end_pos COMMENT
	start_pos=-1
	end_pos=-1
	q_start_pos=-1
	q_end_pos=-1
	line_start=-1
	stop_parsing=0
	COMMENT=""
	subj_seq=""
	query_seq=""
	if ( /^ Score/) (line_start=NR)

}

if ( /^ Score/)
{	
	
	split($0, LINE_SPLIT, "=")
	split(LINE_SPLIT[3], LINE_SPLITTED, ",")
	while (sub(/ /, "", LINE_SPLITTED[1])){}
	evalue = LINE_SPLITTED[1]
	
	if (!/dentities/){getline}
	split($0, LINE_SPLIT, "%")
	
	while (sub(/[\(\)QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm=\:]/, "", LINE_SPLIT[1])){}
	split(LINE_SPLIT[1], LINE_SPLITTED, " ")
	identity = LINE_SPLITTED[length(LINE_SPLITTED)]


}


  if ( /^>/)   {
	subject = substr($1, 2)

	if (line_start==-1)
	{line_start = NR }


}






}
