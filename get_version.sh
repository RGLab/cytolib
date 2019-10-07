sed -n -r '/^Version: /{s/^Version:[ \t]+([^ \t]+)[ \t]*$/\1/;p}' DESCRIPTION
