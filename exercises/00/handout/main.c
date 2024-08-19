#include <stdlib.h>
#include <stdio.h>
#include "bitmap.h"

#define XSIZE 2560 // Size of before image
#define YSIZE 2048



int main(void) {
    uchar *image = calloc(XSIZE * YSIZE * 3, 1); // Three uchars per pixel (RGB)
    if (image == NULL) {
        printf("Error allocating memory for image.\n");
        return 1;
    }

    readbmp("before.bmp", image);

    // Invert the image (upside down)
    invert_image(image, XSIZE, YSIZE);

    // Resize the image
    uchar *resized_image = NULL;
    resize_image(image, XSIZE, YSIZE, &resized_image);

    // Invert the colors of the resized image
    invert_colors(resized_image, XSIZE * 2, YSIZE * 2);

    savebmp("after.bmp", resized_image, XSIZE * 2, YSIZE * 2);

    free(image);
    free(resized_image);
    return 0;
}
