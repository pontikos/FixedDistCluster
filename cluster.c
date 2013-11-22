/*
 * =====================================================================================
 *
 *       Filename:  cluster.c
 *
 *    Description: Space efficient fixed-distance 3D clustering algorithm.  Uses region growing approach.
 *
 *        Version:  1.0
 *        Created:  03/08/09 01:51:07 BST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Nikolas Pontikos, Msci Bioinformatics, MEng Computer Science
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#ifdef DEBUG
#define DBG(x) printf x
#else
#define DBG(x)
#endif

typedef struct voxel {
	int coordinate[3];
	int cluster;
} VOXEL;

typedef struct linkedlist {
	struct linkedlist *next;
	VOXEL *voxel;
} LINKEDLIST;

unsigned int total_number_of_elements = 0;
unsigned int elements_assigned_to_clusters = 0;

LINKEDLIST* read_voxels(char* filename, bool header) {

	FILE *file = fopen ( filename, "r" );

	if ( file == NULL ) perror ( filename );

	char line [ 128 ];
	const char const delim[] = ",";
	const char const newline[] = "\n";
	LINKEDLIST *item = NULL, *previous_item = NULL;

	//skip first line
	if (header) fgets ( line, sizeof line, file );

	while ( fgets ( line, sizeof line, file ) != NULL ) /* read a line */ {
		//strip newline
		//fputs ( line, stdout ); /* write the line */
		line[strcspn(line,newline)] = '\0';
		char* tok = NULL;
		tok = strtok(line, delim);
		if (tok == NULL) {
			continue;
		}
		item = malloc(sizeof(LINKEDLIST));
		item->voxel = malloc(sizeof(VOXEL));
		item->next = previous_item;
		item->voxel->coordinate[2] = atoi(tok);
		tok = strtok(NULL, delim);
		item->voxel->coordinate[1] = atoi(tok);
		tok = strtok(NULL, delim);
		item->voxel->coordinate[0] = atoi(tok);
		item->voxel->cluster = -1;
		total_number_of_elements++;
		previous_item = item;
	}
	fclose ( file );
	return item;
}

float euclidean_distance(const VOXEL const * p, const VOXEL const * p2) {
	if (sizeof(p->coordinate) != sizeof(p2->coordinate)) perror("different size coordinates");
	float dist = 0;
	int i = 0;
	for (i=0; i < sizeof(p->coordinate)/sizeof(int); i++) dist += pow((p->coordinate[i]-p2->coordinate[i]),2);
	return sqrt(dist);
}

#ifdef DEBUG
int count_items(LINKEDLIST const * head, bool print) {
	int c = 0;
	while (head != NULL) {
		c++;
		if (print == true) printf("%d,%d,%d;\n", head->voxel->coordinate[0],head->voxel->coordinate[1],head->voxel->coordinate[2]);
		head = head->next;
	} 
	return c;
}
#endif

//deletes next item and returns this one
LINKEDLIST** delete_next_item(LINKEDLIST **pe) {
	LINKEDLIST *e = (*pe)->next;
	if (e == NULL) return pe;
	(*pe)->next = e->next;
	free(e);
	return pe;
}

//deletes this item and returns next one
LINKEDLIST* delete_head(LINKEDLIST *h) {
	if (h == NULL) return h;
	LINKEDLIST *nh = h->next;
	free(h);
	return nh;
}

LINKEDLIST* get_neighbours(const VOXEL const * p, LINKEDLIST **h, const float distance) {
	LINKEDLIST *head = *h;
	LINKEDLIST *pe = head, *e = head->next;
	LINKEDLIST *neighbour = NULL, *previous_neighbour = NULL;
	DBG(("Getting neighbours of %d,%d,%d\n", p->coordinate[0], p->coordinate[1], p->coordinate[2]));
	while (e != NULL) {
		//this should never be the case
		//if (e->voxel->cluster != -1) exit(0);
		if (euclidean_distance(p, e->voxel) <= distance) {
			DBG(("neighbour %d,%d,%d\n", e->voxel->coordinate[0], e->voxel->coordinate[1], e->voxel->coordinate[2]));
			printf("%u\n", total_number_of_elements-(++elements_assigned_to_clusters));
			e->voxel->cluster = p->cluster;
			//create and append neighbour
			neighbour = malloc(sizeof(LINKEDLIST));
			neighbour->next = previous_neighbour;
			neighbour->voxel = e->voxel;
			previous_neighbour = neighbour;
			//remove from pointstovisit
			pe->next = e->next;
			free(e);
			e = pe->next;
		} else {
			pe = e;
			e = e->next;
		}
	}
	//free head if need be, head could be same as p
	if (head->voxel->cluster == -1 && euclidean_distance(p, head->voxel) <= distance) {
		DBG(("head neighbour %d,%d,%d of %d,%d,%d\n", head->voxel->coordinate[0], head->voxel->coordinate[1], head->voxel->coordinate[2], p->coordinate[0], p->coordinate[1], p->coordinate[2]));
		printf("%u\n", total_number_of_elements-(++elements_assigned_to_clusters));
		head->voxel->cluster = p->cluster;
		//create and append neighbour
		neighbour = malloc(sizeof(LINKEDLIST));
		neighbour->next = previous_neighbour;
		neighbour->voxel = head->voxel;
		//remove from pointstovisit
		head = delete_head(head);
	}
	*h = head;
	return neighbour;
}

int bfcluster(LINKEDLIST* pointstovisit, const float distance) {
	int cluster = -1;
	LINKEDLIST *neighbour=NULL, *new_neighbour=NULL, *tail=NULL;
#ifdef DEBUG
	int c = 0;
#endif
	while (pointstovisit != NULL) {
		pointstovisit->voxel->cluster = ++cluster;
		DBG(("STARTING CLUSTER %d\n", cluster));
		DBG(("v %d,%d,%d;\n", pointstovisit->voxel->coordinate[0],pointstovisit->voxel->coordinate[1],pointstovisit->voxel->coordinate[2]));
		neighbour = get_neighbours(pointstovisit->voxel, &pointstovisit, distance);
		tail = NULL;
		while (neighbour != NULL && pointstovisit != NULL) {
			DBG(("\nNeighbours of %d,%d,%d\n", neighbour->voxel->coordinate[0],neighbour->voxel->coordinate[1],neighbour->voxel->coordinate[2]));
			new_neighbour = get_neighbours(neighbour->voxel, &pointstovisit, distance);
			if (new_neighbour == NULL) {
				DBG(("No new neighbours\n"));
				neighbour = delete_head(neighbour);
				continue;
			}
#ifdef DEBUG
			c = count_items(new_neighbour, true);
			printf("New neighbours count: %d\n", c);
#endif
			//add to tail of neighbour linked list
			//potential speed up here is to use last tail
			//if (tail == NULL)
				tail = neighbour;
			while (tail->next != NULL) tail = tail->next;
			tail->next = new_neighbour;
			DBG(("Neighbour list count: %d\n",count_items(neighbour, false)));
			DBG(("Points left to visit %d\n",count_items(pointstovisit,false)));
			//DBG(("%d\n", count_items(neighbour, false)+count_items(pointstovisit,false)));
			neighbour = delete_head(neighbour);
		}
		DBG(("FINISHED CLUSTER %d\n", cluster));
		pointstovisit = delete_head(pointstovisit);
	}
	return cluster;
}


void print_clusters(LINKEDLIST* e, int maxcluster, char* outfilename) {
	FILE *outfile = fopen (outfilename , "w" );
	LINKEDLIST* head = e;
	int i, n=0, c;
	for (i=-1; i<=maxcluster; i++) {
		e = head;
		c = 0;
		while (e != NULL) {
			if (e->voxel->cluster == i) {
				char line[128];
				sprintf(line, "%d,%d,%d,%d\n", e->voxel->cluster, e->voxel->coordinate[0], e->voxel->coordinate[1], e->voxel->coordinate[2]);
				fputs(line, outfile);
				c++;
				n++;
			}
			e = e->next;
		}
		printf("Cluster %d contains %d elements\n", i, c);
	}
	fclose(outfile);
	printf("Total # of elements %d\n", n);
}



int main(int argc, char** args) {

	if (argc < 2) {
		printf("Usage: %s <filename>\n", args[0]);
		return 1;
	}

	LINKEDLIST *head = read_voxels(args[1], true);
	LINKEDLIST *v = head;
	//create copy of the linked list
	LINKEDLIST *prev_pointstovisit = NULL, *pointstovisit = NULL;
	while (v != NULL) {
		pointstovisit = malloc(sizeof(LINKEDLIST));
		pointstovisit->voxel = v->voxel;
		pointstovisit->next = prev_pointstovisit;
		prev_pointstovisit = pointstovisit;
		v = v->next;
	}
	int maxcluster = bfcluster(pointstovisit, sqrt(2));

	char outfilename[80];
	sprintf(outfilename, "clusters_%s", args[1]);
	printf("Writing clusters to %s\n", outfilename);
	print_clusters(head, maxcluster, outfilename);

	//testing euclidean_distance
	/*--------------------------------------------------
	* VOXEL v;
	* v.coordinate[0] = 1;
	* v.coordinate[1] = 1;
	* v.coordinate[2] = 1;
	* VOXEL v2;
	* v2.coordinate[0] = 0;
	* v2.coordinate[1] = 0;
	* v2.coordinate[2] = 0;
	* printf("%f\n", euclidean_distance(&v2, &v));
	*--------------------------------------------------*/

	return 0;
}


