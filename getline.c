#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>

int mygetchar() {
    int c;
    while  ((c=getchar()) == '\r') {
    }
    return c;
}

int get_line(char oneline[]) {
    int c;
    int len = 0;
    while ((c=mygetchar() != EOF)) {
        if (c == '\n') {
            break;
        }
        oneline[len] = c;
        len++;
    }
    oneline[len] = '\0';
    return 0;
}

int main(int argc, char** argv) {
    char oneline[999];
    get_line(oneline);
    printf("%c\n", oneline[1]);
    return 0;
}