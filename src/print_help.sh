# This script is meant to be called with a string to grep with; however
# If no string is provided, the program should still work
grep_string=$1
echo $grep_string

# Find the length of the helpfile, subtract the license (37 lines)
num_lines=$(cat src/helpfile | wc -l)
license_lines=3
print_lines=$((num_lines-license_lines))

# print help
if [ -n "$grep_string" ]; then
    tail src/helpfile -n $print_lines | grep $grep_string
else
    tail src/helpfile -n $print_lines
fi

