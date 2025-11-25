#include <stdio.h>

int main() 
{
    FILE *parameter_file = fopen("/home/cillian/.eff_source_J/notebooks/parameters.dat", "r");
    char line[256];
    double parameters[5];
    fgets(line, sizeof(line), parameter_file);
    double rmin, rmax, e, l, Trad, halfTrad;
    int arr_length;
    sscanf(line, "%lf %lf %lf %lf %d %lf %lf", &rmin, &rmax, &e, &l, &arr_length, &Trad, &halfTrad);
    fclose(parameter_file);
    printf("%d", arr_length);
    return 0;
}
