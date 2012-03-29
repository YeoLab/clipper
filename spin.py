import sys
import time
import os

def spinning_cursor():
    cursor='/-\|'
    i = 0
    while 1:
        yield cursor[i]
        i = (i + 1) % len(cursor)
def spin(filename):
    print "waiting for file: %s to be created" %(filename)
    for c in spinning_cursor():
        sys.stdout.write(c)
        sys.stdout.flush()
        time.sleep(10)
        sys.stdout.write('\b')
        if os.path.exists(filename):
            break
    print "All done"


if __name__ == "__main__":
    spin(sys.argv[1])
