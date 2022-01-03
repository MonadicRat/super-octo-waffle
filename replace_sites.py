import re

meth = open("/media/yanus/3ce7ab71-4bf7-4250-894c-508fb7a28109/Meth/coeffs.sample.txt")
for line in meth:
    print(re.sub(r' {1,}', '\t', line), end='')
