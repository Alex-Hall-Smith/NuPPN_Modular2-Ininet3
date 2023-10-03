#!/usr/bin/gawk -f
/^#(ifdef|ifndef)\>/ {matches[$2, FILENAME] = ""}
/^#if\>/ {
        founderr = 0 # set to 1 when an error was found
        ind = 0
        while(1) {
                match(substr($0, ind), /defined\(([^)]+)\)/, m)
                if (RSTART == 0) break
                matches[m[1], FILENAME] = ""
                ind += RSTART + RLENGTH
        }
}
{
        if (!match(FILENAME, /\.[fF]90$/)) next
        line = gensub(/\\\\|\\'|\\"/, "", "g") # remove escape chars
        gsub(/'.*'/, "", line) # remove strings
        gsub(/".*"/, "", line) # remove strings
        sub( /!.*/,  "", line) # remove comments
        gsub( /\.[[:alpha:]]+\./, "", line) # remove fortran operators
        if (     match(line, /[[:digit:]]*\.[[:digit:]]+(e[+-]?[[:digit:]]+)?\>/, m) \
              || match(line, /[[:digit:]]+\.[[:digit:]]*(e[+-]?[[:digit:]]+)?\>/, m) \
              || match(line, /[[:digit:]]+e[+-]?[[:digit:]]+\>/, m) \
           )
        {
                start = index($0,m[0])
                s = sprintf("%s:%d: ", FILENAME, FNR)
                print s $0
                for (i=0; i<length(s)-1+start; i++) printf " "
                printf "^"
                for (i=0; i<RLENGTH-2; i++) printf " "
                print "|"
                founderr = 1
        }
}
END {
        while (1) {
                gotline = getline < "checker-whitelist"
                if (gotline <= 0) break
		if ($0 !~/^#/) whitelist[$0] = ""
        }
        for (i in matches) {
                split(i, a, SUBSEP)
                if (!(a[1] in whitelist)) {
                        print a[2] ": " a[1]
                        founderr = 1
                }
        }

        exit founderr
}
