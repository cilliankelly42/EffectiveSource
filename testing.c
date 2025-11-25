#include <stdio.h>
#include <sys/stat.h>

int main()
{
    double e = 0.7;
    char path[256];
    sprintf(path, "./eccentric_%.2f", e);
    mkdir(path, 0777);
    return 0;
}


