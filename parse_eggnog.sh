#! /bin/bash
#
# Peter W. Schafran ps997@cornell.edu
#
# Usage: parse_eggnog.sh QUERY_SEQ_NAME[.txt] EGGNOG_COLUMN1 EGGNOG_COLUMN2...EGGNOG_COLUMNn
#
# Accepted EGGNOG_COLUMNS:
# SEED_ORTHOLOG, SEED_EVALUE, SEED_SCORE, BEST_TAX_LEVEL, PREFERRED_NAME, GO,
# EC, KEGG_KO, KEGG_PATHWAY, KEGG_MODULE, KEGG_REACTION, KEGG_RCLASS, BRITE,
# KEGG_TC, CAZY, BIGG, TAXONOMY, EGGNOG_OG, BEST_EGGNOG_OG, COG, FREE_TEXT
# or use ALL for whole line

#help=```
#For eggnog-mapper v2.1.5 to 2.1.10.
#
#Usage: parse_eggnog.sh query_sequence_id[.txt] [prefix].emapper.annotations column1 column2...columnN
#
#query_sequence_id can be a string of the name of a sequence, or a text file with one sequence name per line.
#
#prefix.emapper.annotations is the eggnog-mapper output file off that name.
#
#You can specify as many columns as you want from the headings below, or "ALL" for the whole line.
#
#Accepted column names: seed_ortholog, evalue, score, eggnog_ogs, max_annot_lvl, cog_category, description, preferred_name, gos, ec, kegg_ko, kegg_pathway, kegg_module, kegg_reaction, kegg_rclass, brite, kegg_tc, cazy, bigg_reaction, pfams
#```



QUERY=$1
ANNOTATIONS=$2
PARAMS=()

# Check for help flag
if $1=="-h"
	then echo $help
	exit(0)
elif $1=="--help"
	then echo $help
	exit(0)
fi

# Determine input is file or string
if [ -f $1 ]
  then GREP_TYPE="-f"
else
  GREP_TYPE="-w"
fi

if [ ! -f $2 ]
	then echo "Eggnog annotation file not found"
	exit(1)
fi

if [ $3 = "ALL" ]
  then grep -F $GREP_TYPE "$QUERY" $ANNOTATIONS
else
  # Convert text command line args to corresponding column number in eggnog file
  for i in "${@:3}"
    do if [ $i = "seed_ortholog" ]
      then PARAMS+=(2)
    elif [ $i = "evalue" ]
      then PARAMS+=(3)
    elif [ $i = "score" ]
      then PARAMS+=(4)
    elif [ $i = "eggnog_ogs" ]
      then PARAMS+=(5)
    elif [ $i = "max_annot_lvl" ]
      then PARAMS+=(6)
    elif [ $i = "cog_category" ]
      then PARAMS+=(7)
    elif [ $i = "description" ]
      then PARAMS+=(8)
    elif [ $i = "preferred_name" ]
      then PARAMS+=(9)
    elif [ $i = "gos" ]
      then PARAMS+=(10)
    elif [ $i = "ec" ]
      then PARAMS+=(11)
    elif [ $i = "kegg_ko" ]
      then PARAMS+=(12)
    elif [ $i = "kegg_pathway" ]
      then PARAMS+=(13)
    elif [ $i = "kegg_module" ]
      then PARAMS+=(14)
    elif [ $i = "kegg_reaction" ]
      then PARAMS+=(15)
    elif [ $i = "kegg_rclass" ]
      then PARAMS+=(16)
    elif [ $i = "brite" ]
      then PARAMS+=(17)
    elif [ $i = "kegg_tc" ]
      then PARAMS+=(18)
    elif [ $i = "cazy" ]
      then PARAMS+=(19)
    elif [ $i = "bigg_reaction" ]
      then PARAMS+=(20)
    elif [ $i = "pfams" ]
      then PARAMS+=(21)
    fi
  done

  # Remove any double entries in input
  UNIQ_PARAMS=($(printf "%s\n" "${PARAMS[@]}" | sort | uniq))

  # Convert array to comma separated list to pass to cut
  DELIM=""
  FIELDS=""
  for i in "${UNIQ_PARAMS[@]}"
    do FIELDS="$FIELDS$DELIM$i"
    DELIM=","
  done

  grep -F $GREP_TYPE "$QUERY" $ANNOTATIONS | cut -f 1,"$FIELDS"
fi
