#include "wsp.h"
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
/************/
/*** Main ***/
/************/

std::ofstream output_file;

int main(int argc, char *argv[]){
  int rows{-1}; // counting rows, must be less than num to start
  int num{0}, dim, i, j, num_wsp;
  double sep_const;
  tree_node *root;
  point *pts;
  list_set orders;

  Timer timer;
  if (argc < 3) {
    fprintf(stdout, "Usage 1: wsp num_points num_dim sep_constant\n");
    fprintf(stdout, "Usage 2: wsp file_name sep_constant\n");

    exit(-1);
  }

  if (argc == 4)
  {
    num = atoi(argv[1]);
    dim = atoi(argv[2]);
    sep_const = atof(argv[3]);
    /* generate points */
    /*  printf("Generating points...\n"); */
    /* pts = GeneratePointsDave(num, dim); */
    printf("Genereting points... ");
    timer.start();
    pts = GeneratePoints(num, dim);
    timer.stop();
    cout << " done in " << timer.msecs() << " msecs" << endl;

  }

  string file_name ;
  
  if (argc == 3) // read from the file 
  {
    file_name = argv[1];
    string output_file_name = string(file_name) + ".wsp_out.txt";
    cout << "wsp: " << file_name << " --> " << output_file_name << endl;
    //output_file.open(file_name + ".output.txt");
    output_file.open(output_file_name);
    sep_const = atof(argv[2]);
    // read file to pts data structure

    cout << "Reading file: " << file_name << ", target sep: " << sep_const << endl;
    timer.start();
    ifstream file(file_name);
    bool first_line = true;
    vector<double> coords;
    
    //int rows  = 0; // counting rows
    if (file.is_open()) {
      // EDIT:
      while (!file.eof() && rows < num) {
        string line;
        getline (file, line); // reads a line
        //cout << line << endl;
        istringstream iss(line);
        // EDIT: peel of the header comment marker
        // TODO: Make sure the file is read right
        if (iss.peek() == '#') iss.get();
        if (first_line)
        {
          rows = 0;
          iss >> num;
          iss >> dim;
          first_line = false;
          
          // allocate points data-structure!!!
          pts = (point *) malloc(sizeof(point) * num);	
          for (i = 0; i < num; i++) {
            pts[i].coord = (double *) malloc(sizeof(double)*dim);
            pts[i].lists = (list_elt **) malloc(sizeof(list_elt *) * dim);
          }
          cout << "nr_points = " << num << " and dim = " << dim << endl;
        }else
        {
          //cout << "first_line = false" << endl;
          for (j=0; j<dim; j++) 
          {
            double  sub;
            iss >> sub; // parting the each number in the line
            pts[rows].coord[j] = sub;
            //cout << "part-string: (i,j): " <<rows << ", " << j << " : " <<  sub << endl;
          }
          pts[rows].index = rows;
          rows++;
        }
      }
      file.close();
    }

    timer.stop();
    cout << "done in " << timer.secs() << " secs" <<  endl;          
    
  }
    
  
  orders = (list *)malloc(sizeof(list) * dim);


  /* generate sorted lists */
  /* printf("Sorting points...\n"); */
  printf("Sorting points ... ");
  timer.start();
  for(i = 0; i < dim; i++)
    SortPoints(pts, num, i, &(orders[i]));
  timer.stop();
  cout << "done in " << timer.secs() << " secs" <<  endl;
  
  /* recursive call to BuildTree */
  /* printf("Building fair split-tree...\n"); */
  printf("Building Tree ... ");
  timer.start();
  root = BuildTree(orders, num, dim);
  timer.stop();
  cout << "done in " << timer.secs() << " secs" <<  endl;
  
  /* PrintTree(root, dim, 0); */

  /* recursive call to FindWSP */
  /* printf("Finding well-separated pairs...\n"); */
  cout << "Finding pairs ..." << flush;
  timer.start();
  num_wsp = FindWSP(root, sep_const, dim);
  timer.stop();
  cout << "done in " << timer.secs() << " secs" <<  endl;  
  printf("%lf, %d;\n", sep_const, num_wsp);
  cout << "Upper bound on the number of WS pairs: " << upper_bound(num, dim, sep_const) << endl;
}

/*****************/
/*** FindSplit ***/
/*****************/

void *FindSplit(list_set ls, int dim, 
		half *smaller_half, int *max_dim, int *num_smaller) {

  int i;
  double lmax, lmax_half;
  bool found;
  list_elt *lower_elt_ptr, *upper_elt_ptr;
  point p1, p2;
  
  /* finding the maximum dimension */
  lmax = -1;
  for (i=0; i<dim; i++) {
    p1 = *(ls[i].first->pt);
    p2 = *(ls[i].last->pt);
    if (p2.coord[i] - p1.coord[i] > lmax) {
      lmax = p2.coord[i] - p1.coord[i];
      *max_dim = i;
    }
  }
  p1 = *(ls[*max_dim].first->pt);
  lmax_half = p1.coord[*max_dim] + lmax/2.0;
  
  /* finding the midpoint in dimension max_dim */
  found = false;
  *num_smaller = 0;
  lower_elt_ptr = ls[*max_dim].first;
  upper_elt_ptr = ls[*max_dim].last;
  p1 = *(lower_elt_ptr->pt);
  p2 = *(upper_elt_ptr->pt);

  if (p1.coord[*max_dim] > lmax_half) {
    found = true;
    *smaller_half = lower;
  }
  if ((p2.coord[*max_dim] < lmax_half) && (found == false)) {
    found = true;
    *smaller_half = upper;
  }

  while (found == false) {
    (*num_smaller)++;
    lower_elt_ptr = lower_elt_ptr->next;
    upper_elt_ptr = upper_elt_ptr->prev;
    p1 = *(lower_elt_ptr->pt);
    p2 = *(upper_elt_ptr->pt);

    if (p1.coord[*max_dim] > lmax_half) {
      found = true;
      *smaller_half = lower;
    }
    if ((p2.coord[*max_dim] < lmax_half) && (found == false)) {
      found = true;
      *smaller_half = upper;
    }
  }
}

/*******************/
/*** DeletePoint ***/
/*******************/
void DeletePoint(int dim, point *p, list_set ls) {
  int i;
 
  /* using the set of dim lists in the point to delete the point */
  /* from all dim lists */

  for (i=0; i<dim; i++) {
    if(p->lists[i]->prev != NULL)
      p->lists[i]->prev->next = p->lists[i]->next;
    else {
      ls[i].first = p->lists[i]->next;
	    p->lists[i]->next->prev = NULL;
    }
    if(p->lists[i]->next != NULL)
      p->lists[i]->next->prev = p->lists[i]->prev;
    else {
      ls[i].last = p->lists[i]->prev;
      p->lists[i]->prev->next = NULL;
    }
  }
}

/******************/
/*** FindCenter ***/
/******************/
void FindCenter(list_set ls, tree_node *tnode, int dim) {
  int i, max_dim;
  double lmax, half, distsqrd = 0;
  point p1, p2;
  
  lmax = -1;
  for (i=0; i<dim; i++) {
    p1 = *(ls[i].first->pt);
    p2 = *(ls[i].last->pt);
    tnode->center[i] = (p2.coord[i] + p1.coord[i]) / 2;
    half = ((p2.coord[i] - p1.coord[i])/2.0);
    distsqrd += (half*half);
    if (p2.coord[i] - p1.coord[i] > lmax) {
      lmax = p2.coord[i] - p1.coord[i];
      max_dim = i;
    }
  }
  
  tnode->radius = sqrt(distsqrd);
}

void FreeListSet(list_set ls, int dim) {
  int i;
  
  for(i = 0; i < dim; i++) 
    free(ls[i].mem);
  free(ls);
}

/*****************/
/*** BuildTree ***/
/*****************/

tree_node *BuildTree(list_set ls, int num, int dim){
  int i, j, total;
  list_set ls_copy;

  tree_node *current, *tnode;
    
  int max_dim, num_smaller, num_new_lists;
  half smaller_half;
  point *pp1;
  list_set_elt *new_lists_first, *new_lists_current, *to_free;
  list_elt *new_list_elts, *list_elt_current;

  /* base case */
  if (num == 1) {
    tnode = (tree_node *)malloc(sizeof(tree_node));
    tnode->center = (double *)malloc(sizeof(double) * dim);
    tnode->lchild = NULL;
    tnode->rchild = NULL;
    tnode->pt = ls[0].first->pt;
    tnode->nr_pt = 1;
    for(i = 0; i < dim; i++) 
      tnode->center[i] = ls[0].first->pt->coord[i];
    tnode->radius = 0;
    //for(i = 0; i < dim; i++) 
      //free(ls->first);
    free(ls);
    return (tnode);
  }

  /* copying the lists */
  ls_copy = (list *)malloc(sizeof(list) * dim);
  for (i=0; i<dim; i++) {
    ListCopy(&(ls[i]), &(ls_copy[i]), num);
    ls[i].mem = ls[i].first;
    ls_copy[i].mem = ls_copy[i].first;
  }

  num_new_lists = 0; 
  
  /* main loop */
  total = num;
  while (total > num/2) {

    /* create tree_node for this one, and fill in radius, center */
    if(num_new_lists == 0) {
      tnode = current = (tree_node *)malloc(sizeof(tree_node));
    }
    else {
      current->lchild = (tree_node *)malloc(sizeof(tree_node));
      current = current->lchild;
    } 
    current->pt = NULL;
    current->lchild = current->rchild = NULL;
    current->center = (double *)malloc(sizeof(double) * dim);
    current->nr_pt = total; // set the number of points

    FindCenter(ls, current, dim);

    /* find the widest dimension and split along it */
    FindSplit(ls, dim, &smaller_half, &max_dim, &num_smaller);

    /* creating the duplicate lists */
    num_new_lists++;
    if(num_new_lists == 1){ 
      new_lists_first = new_lists_current 
	= (list_set_elt *) malloc(sizeof(list_set_elt));
      new_lists_current->next = NULL;
    }
    else {
      new_lists_current->next = (list_set_elt *)malloc(sizeof(list_set_elt));
      new_lists_current = new_lists_current->next;
      new_lists_current->next = NULL;
    }
    
    new_lists_current->data = (list *)malloc(sizeof(list)*dim);
    new_lists_current->count = num_smaller;

    for (i=0; i<dim; i++) {
      new_list_elts = (list_elt *) malloc(sizeof(list_elt)*num_smaller);
      DoubleLink(new_list_elts, num_smaller, &(new_lists_current->data[i]));
    }
    
    total -= num_smaller;

    /* deleting elements in the smaller half from all existing lists */
    if (smaller_half == lower) {
      list_elt_current = ls[max_dim].first;
      while (num_smaller != 0) {
	pp1 = list_elt_current->pt;

	/* setting the pointer to the new set of lists that the point */
	/* will be added to */
	pp1->new_list = new_lists_current->data;

	DeletePoint(dim, pp1, ls);

	num_smaller--;
	list_elt_current = list_elt_current->next;
      }
    }
    if (smaller_half == upper) {
      list_elt_current = ls[max_dim].last;
      while (num_smaller != 0) {
	pp1 = list_elt_current->pt;

	pp1->new_list = new_lists_current->data;
	
	DeletePoint(dim, pp1, ls);
	
	num_smaller--;
	list_elt_current = list_elt_current->prev;
      }
    } 
  }

  /* create a copy for the remaining points, so we can free memory */
  num_new_lists++;
  if(num_new_lists == 1){ 
    new_lists_first = new_lists_current 
      = (list_set_elt *) malloc(sizeof(list_set_elt));
    new_lists_current->next = NULL;
  }
  else {
    new_lists_current->next = (list_set_elt *)malloc(sizeof(list_set_elt));
    new_lists_current = new_lists_current->next;
    new_lists_current->next = NULL;
  }
  new_lists_current->data = (list *)malloc(sizeof(list) * dim);
  new_lists_current->count = total;
  for (i=0; i<dim; i++) {
    ListCopy(&(ls[i]), &(new_lists_current->data[i]), total);
    list_elt_current = new_lists_current->data[i].first;
    for(j = 0; j<total; j++) {
      list_elt_current->pt->lists[i] = list_elt_current;
      list_elt_current = list_elt_current->next;
    }
  }

  /* now we can free the original list */
  FreeListSet(ls, dim);

  new_lists_current = new_lists_first;
  for(i = 0; i < num_new_lists; i++) {
    for(j = 0; j < dim; j++) 
      new_lists_current->data[j].last = new_lists_current->data[j].first;
    new_lists_current = new_lists_current->next;
  }

  /* run through copy of original lists, for each point, insert into */
  /* the correct new list */
  for(j = 0; j < dim; j++) {
    list_elt_current = ls_copy[j].first;
      for(i = 0; i < num; i++) {
	if(list_elt_current->pt->new_list != NULL) {
	  list_elt_current->pt->new_list[j].last->pt = list_elt_current->pt;
	  list_elt_current->pt->lists[j] 
	    = list_elt_current->pt->new_list[j].last;
	  list_elt_current->pt->new_list[j].last
	    = list_elt_current->pt->new_list[j].last->next;
	}

	list_elt_current = list_elt_current->next;
      }
  }

  /* now we can free the copies of the lists */
  FreeListSet(ls_copy, dim);

  /* restore "last" pointers for each list */
  new_lists_current = new_lists_first;
  for(i = 0; i < num_new_lists; i++) {
    for(j = 0; j < dim; j++) {
      list_elt_current = new_lists_current->data[j].first;
      while(list_elt_current->next)
	list_elt_current = list_elt_current->next;
      new_lists_current->data[j].last = list_elt_current;
    }
    new_lists_current = new_lists_current->next;
  }

  /* recursive call */
  current = tnode;
  new_lists_current = new_lists_first;
  for(i = 0; i < num_new_lists - 1; i++) {
    current->rchild = 
      BuildTree(new_lists_current->data, new_lists_current->count, dim);
    if(i != (num_new_lists - 2)) 
      current = current->lchild;
    else
      current->lchild = BuildTree(new_lists_current->next->data, 
				  new_lists_current->next->count, dim);
    to_free = new_lists_current;
    new_lists_current = new_lists_current->next;
    free(to_free);
  }

  /* return tree node */
  return (tnode);
}
