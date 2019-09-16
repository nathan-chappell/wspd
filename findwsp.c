#include <fstream>
#include <assert.h>
#include <vector>
#include "wsp.h"

int FindWSP2(tree_node *tnode1, tree_node *tnode2, double s, int dim);
int wellsep(tree_node *tnode1, tree_node *tnode2, double s, int dim, double& dist);

/***************/
/*** FindWSP ***/
/***************/

int FindWSP(tree_node *tnode, double s, int dim){

  if (tnode->pt != NULL) {
    /* printf("tnode->pt != NULL!!!\n"); */
    return (0);
  }
  else {
    return (FindWSP(tnode->lchild, s, dim)
	    + FindWSP(tnode->rchild, s, dim)
	    + FindWSP2(tnode->lchild, tnode->rchild, s, dim)); 
  }
}

int dfs(tree_node *tnode, vector<double*>& pt)
{
	if (tnode->lchild == NULL && tnode->rchild == NULL){
		pt.push_back(tnode->pt->coord);
		return 1;
	}
	return dfs(tnode->lchild, pt) + dfs(tnode->rchild, pt);
}

// EDIT: get index information
using Points = vector<point*>;

int dfs(tree_node *tnode, Points& pts)
{
	if (tnode->lchild == NULL && tnode->rchild == NULL){
		pts.push_back(tnode->pt);
		return 1;
	}
	return dfs(tnode->lchild, pts) + dfs(tnode->rchild, pts);
}


/****************/
/*** FindWSP2 ***/
/****************/

// EDIT:
// TODO: output the information about the nodes
//
// idea: one file has: indices of pairs
//       other file has: separation, radius

// format: [[left indices],[right indices]]
std::ofstream dumbell_i_file{"dumbells.txt"};

void output_indices(Points &l, Points &r) {
  for (auto p : l) dumbell_i_file << p->index << " ";
  dumbell_i_file << " | ";
  for (auto p : r) dumbell_i_file << p->index << " ";
  dumbell_i_file << "\n";
}

// format: sep l_n, r_n, l_c, r_c, l_r, r_r
std::ofstream dumbell_s_c_r_file{"dumbells_s_c_r.txt"};

//TODO: make sure we don't deref a nullptr
void output_info(double sep, tree_node *l, tree_node *r) {
  dumbell_s_c_r_file << sep << " "
            //<< l->nr_pt << " "
            //<< r->nr_pt << " "
            << *l->center << " "
            << *r->center << " "
            << l->radius << " "
            << r->radius << "\n";
}

int FindWSP2(tree_node *tnode1, tree_node *tnode2, double s, int dim){

  /*   printf("calling FindWSP2!\n"); */
  double distance;
  if (wellsep(tnode1, tnode2, s, dim, distance)) {
	  // file format:
	  // tnode1->nr_pt,tnode2->nr_pt,tnode1->center,tnode2->center,tnode1->radius,tnode2->radius
	  //
    // EDIT:
	  Points tnode1_pts, tnode2_pts;
	  dfs(tnode1, tnode1_pts);
	  dfs(tnode2, tnode2_pts);
    output_indices(tnode1_pts,tnode2_pts);
    output_info(s,tnode1,tnode2);

    /*
	  vector<double*> tnode1_pts, tnode2_pts;
	  dfs(tnode1, tnode1_pts);
	  dfs(tnode2, tnode2_pts);

	  assert(tnode1_pts.size() == tnode1->nr_pt);
	  assert(tnode2_pts.size() == tnode2->nr_pt);

	  std::ofstream outfile;
	  outfile.open("wsp.out", std::ios_base::app);
	  outfile<< tnode1->nr_pt << "," << tnode2->nr_pt << "," << distance << ",";
	  for (int j=0;j<tnode1->nr_pt; j++)
		  for (int i=0;i<dim;i++)
        	          outfile << tnode1_pts[j][i] << ",";
          for (int j=0;j<tnode2->nr_pt; j++)
	        for (int i=0;i<dim;i++)
		{
        		outfile << tnode2_pts[j][i] ;
                        if (i<(dim-1))
                                outfile << ",";
		}

	  outfile << endl;
    */
	 /* outfile<< tnode1->nr_pt << "," << tnode2->nr_pt << ",";
	  for (int i=0;i<dim;i++)
		  outfile << tnode1->center[i] <<",";
	  for (int i=0;i<dim;i++)
		outfile << tnode2->center[i] << ",";
	  outfile << tnode1->radius << ",";
	  outfile << tnode2->radius << endl;
	  */
	  //outfile.close();
	  return (1);
  }
  else {
    if (tnode1->radius > tnode2->radius)
      return (FindWSP2(tnode1->lchild, tnode2, s, dim) 
	      + FindWSP2(tnode1->rchild, tnode2, s, dim));
    else
      return (FindWSP2(tnode1, tnode2->lchild, s, dim) 
	      + FindWSP2(tnode1, tnode2->rchild, s, dim));
  }
}

/***************/
/*** wellsep ***/
/***************/

int wellsep(tree_node *tnode1, tree_node *tnode2, double s, int dim, double& dist){
  int i;
  double radius, distance;

  radius = std::max(tnode1->radius, tnode2->radius);

  distance = 0.0;
  for (i=0; i<dim; i++) {
    distance += (tnode1->center[i] - tnode2->center[i])*
      (tnode1->center[i] - tnode2->center[i]);
  }

  dist = distance = sqrt(distance);
  distance = distance - 2*radius;

  if (distance > s*radius)
    return (1);
  else
    return (0);
}
