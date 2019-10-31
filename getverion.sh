sed -n -e "s/^Version:[ '\\t']\\{1,\\}\\([^ '\\t']\\{1,\\}\\)[ '\\t']*$/\\1/p" DESCRIPTION
