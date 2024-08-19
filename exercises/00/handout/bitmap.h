#ifndef BITMAP_H
#define BITMAP_H


typedef unsigned char uchar;

void savebmp(char *name, uchar *buffer, int x, int y);
void readbmp(char *filename, uchar *array);
void invert_colors(uchar *image, int width, int height);  // New function declaration
void resize_image(uchar *image, int width, int height, uchar **new_image);  // New function declaration
void invert_image(uchar *image, int width, int height);  // New function declaration

#endif
