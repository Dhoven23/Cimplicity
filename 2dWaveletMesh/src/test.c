#include "header.h"

int main()
{
	int i = 0;
	struct timespec ts;
	ts.tv_sec = 0;
	ts.tv_nsec = 10000000;
	while(i < 1000)
	{
		i++;
		printf("%3i", i);


		nanosleep(&ts,&ts);
		printf("\b\b\b");
		printf("%3i", i);
	}
	printf("\n");
	return 0;
}