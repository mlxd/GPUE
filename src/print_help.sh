################################################################################
# print_help - GPUE: Split Operator based GPU solver for Nonlinear 
# Schrodinger Equation, Copyright (C) 2011-2016, Lee J. O'Riordan 
# <loriordan@gmail.com>, James R. Schloss <jrs.schloss@gmail.com>, 
# Tadhg Morgan, Neil Crowley. All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are 
# met:
# 
# 1. Redistributions of source code must retain the above copyright 
# notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright 
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its 
# contributors may be used to endorse or promote products derived from 
# this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
################################################################################

# This script is meant to be called with a string to grep with; however
# If no string is provided, the program should still work
grep_string=$1
echo $grep_string

# Find the length of the helpfile, subtract the license (37 lines)
num_lines=$(cat src/helpfile | wc -l)
license_lines=37
print_lines=$((num_lines-license_lines))

# print help
if [ -n "$grep_string" ]; then
    tail src/helpfile -n $print_lines | grep $grep_string
else
    tail src/helpfile -n $print_lines
fi

