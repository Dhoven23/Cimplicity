#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define LENGTH 5*5
#define BUFF_LENGTH 9*9

bool Index(int x, int y, int* ind);
void PrintBuffer(char* buff);

int main(){
	char* buff = malloc(BUFF_LENGTH*sizeof(*buff));
	memset(buff,' ',BUFF_LENGTH*sizeof(*buff));
	unsigned count = 5;
	int ind = 0;
	for (int i = 0; i < count; ++i){
		for (int j = 0; j < count; ++j){
			int x = 2*j;
			int y = 2*i;
			if (Index(x,y-1,&ind)){
				buff[ind] = '|';
			}
			if (Index(x,y+1,&ind)){
				buff[ind] = '|';
			}
			if (Index(x-1,y,&ind)){
				buff[ind] = '-';
			}
			if (Index(x+1,y,&ind)){
				buff[ind] = '-';
			}
		}
	}

	PrintBuffer(buff);
	free(buff);
}

bool Index(int x, int y, int* ind){
	int p_ind = x + (9*y);
	if ((p_ind < 81) && (p_ind > -1) && (x > -1) && (y > -1) && (y < 9) && (x < 9)){
		*ind = p_ind;
		return true;
	}
	return false;
}

void PrintBuffer(char* buff){
	unsigned count = 9;
	for (int i = 0; i < count; ++i){
		for (int j = 0; j < count; ++j){
			int ind = 0;
			Index(j,i,&ind);
			printf(" %c ",buff[ind]);
		} printf("\n");
	}
}