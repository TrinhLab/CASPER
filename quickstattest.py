"""Super quick stats of file."""


f = open("/Users/brianmendoza/Dropbox/CASPER/OFF_QUERY.txt")
no_lines = 0

for line in f:
    no_lines += 1
f.close()

print(no_lines)