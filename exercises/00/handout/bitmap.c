#include <stdlib.h>
#include <stdio.h>
#include "bitmap.h"






// save 24-bits bmp file, buffer must be in bmp format: upside-down
void savebmp(char *name, uchar *buffer, int x, int y) {
	FILE *f=fopen(name,"wb");
	if(!f) {
		printf("Error writing image to disk.\n");
		return;
	}
	unsigned int size=x*y*3+54;
	uchar header[54]={'B','M',size&255,(size>>8)&255,(size>>16)&255,size>>24,0,
                    0,0,0,54,0,0,0,40,0,0,0,x&255,x>>8,0,0,y&255,y>>8,0,0,1,0,24,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	fwrite(header,1,54,f);
	fwrite(buffer,1,x*y*3,f);
	fclose(f);
}


// read bmp file and store image in contiguous array
void readbmp(char* filename, uchar* array) {
	FILE* img = fopen(filename, "rb");   //read the file
	uchar header[54];
	size_t result = fread(header, sizeof(uchar), 54, img); // read the 54-byte header

	if (result != 54)
	{
		perror("Error reading file.\n");
		return;
	}

  // extract image height and width from header
	int width = *(int*)&header[18];
	int height = *(int*)&header[22];
	int padding=0;
	while ((width*3+padding) % 4!=0) padding++;

	int widthnew=width*3+padding;
	uchar* data = calloc(widthnew, sizeof(uchar));

	for (int i=0; i<height; i++ ) {
		result = fread( data, sizeof(uchar), widthnew, img);

		if (result != (size_t) widthnew)
		{
			perror("Error reading file.\n");
		}

		for (int j=0; j<width*3; j+=3) {
			array[3 * i * width + j + 0] = data[j+0];
			array[3 * i * width + j + 1] = data[j+1];
			array[3 * i * width + j + 2] = data[j+2];
		}
	}
	fclose(img); //close the file
}

// New function that inverts the colors of the image
void invert_colors(uchar *image, int width, int height) {
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            uchar *pixel = &image[(y * width + x) * 3];
            // Invert the colors
            pixel[0] = 255 - pixel[0];
            pixel[1] = 255 - pixel[1];
            pixel[2] = 255 - pixel[2];
        }
    }
}

// Function to resize the image to double its size
void resize_image(uchar *image, int width, int height, uchar **new_image) {
    int new_width = width * 2;
    int new_height = height * 2;

    // Allocate memory for the new image
    *new_image = (uchar *)calloc(new_width * new_height * 3, 1);
    if (*new_image == NULL) {
        printf("Error allocating memory for resized image.\n");
        return;
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // Get the original pixel
            uchar *pixel = &image[(y * width + x) * 3];

            // Map this pixel to 4 pixels in the new image
            for (int dy = 0; dy < 2; dy++) {
                for (int dx = 0; dx < 2; dx++) {
                    int new_x = x * 2 + dx;
                    int new_y = y * 2 + dy;
                    uchar *new_pixel = &(*new_image)[(new_y * new_width + new_x) * 3];
                    new_pixel[0] = pixel[0];
                    new_pixel[1] = pixel[1];
                    new_pixel[2] = pixel[2];
                }
            }
        }
    }
}

// Function to invert the image (upside down)
void invert_image(uchar *image, int width, int height) { 
    // Allocate memory for the new image
	uchar *new_image = (uchar *)calloc(width * height * 3, 1);
	if (new_image == NULL) {
		printf("Error allocating memory for inverted image.\n");
		return;
	}

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			// Get the original pixel
			uchar *pixel = &image[(y * width + x) * 3];
			// Get the new pixel
			uchar *new_pixel = &new_image[((height - y - 1) * width + x) * 3];
			new_pixel[0] = pixel[0];
			new_pixel[1] = pixel[1];
			new_pixel[2] = pixel[2];
		}
	}

	// Copy the new image to the original image
	for (int i = 0; i < width * height * 3; i++) {
		image[i] = new_image[i];
	}

	// Free the memory
	free(new_image);
}

