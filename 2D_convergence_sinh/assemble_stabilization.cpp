//Pressure jump stabilisation.
#if USE_STAB  
 	std::vector<unsigned int> stab_dofs_cols2;
	std::vector<Real> stab_dofs_vals2;
    for (unsigned int s=0; s<elem->n_sides(); s++)
    {
      if (elem->neighbor(s) == NULL)
      {   
				//Only do something on the interior edges.
      } //if (elem->neighbor(s) == NULL)
			else{

      const Elem* neighbor = elem->neighbor(s);
   		AutoPtr<Elem> side (elem->build_side(s));


      std::vector<Real> node_co_ords_x;
      std::vector<Real> node_co_ords_y;
      for (unsigned int n=0; n<elem->n_nodes(); n++)
      {
          Node *e_node = elem->get_node(n);
          for (unsigned int m=0; m<neighbor->n_nodes(); m++)
          {
            Node *n_node = neighbor->get_node(m);
            if( ( (*n_node)(0)==(*e_node)(0) ) && ( (*n_node)(1)==(*e_node)(1) ) ){
              node_co_ords_x.push_back((*n_node)(0));
              node_co_ords_y.push_back((*n_node)(1));
              }
          }
      }
     Real side_length=pow ( pow((node_co_ords_x[0]-node_co_ords_x[1]),2) + pow((node_co_ords_y[0]-node_co_ords_y[1]),2),0.5) ;

	  std::vector<unsigned int> neighbor_dof_indices_p;
      dof_map.dof_indices(neighbor, neighbor_dof_indices_p,p_var);
      const unsigned int n_neighbor_dofs_p = neighbor_dof_indices_p.size();
	  Real hmax=(*elem).hmax();
      Real hmin=(*elem).hmin();
			//Real vol=(*elem).volume();

	  Real delta=dt*DELTA;
      #if !THREED
	  Real factor=-delta*(hmax)*(hmax);;
      //Real factor=-delta*(side_length);
      #endif
      
      #if THREED
      Real factor=-delta*(hmax*hmax*hmax);
      #endif

      //perf_log.push("push back");
      stab_dofs_cols2.push_back(dof_indices_p[0]);
      stab_dofs_vals2.push_back(-factor);
      stab_dofs_cols2.push_back(neighbor_dof_indices_p[0]);
      stab_dofs_vals2.push_back(factor);
      //perf_log.pop("push back");
			}
		}

	//Lots of mallocs happen during the first assemble due to the unexpected sparsity pattern;
	//The implicit_neighbor_dof does not seem to work ?
	//perf_log.push("kstab");
  DenseMatrix<Number> Kstab2;
  Kstab2.resize(1, stab_dofs_vals2.size());
 	for (int i=0; i < stab_dofs_vals2.size(); i++) {
	 	Kstab2(0,i)=stab_dofs_vals2[i];
	}
 	std::vector<unsigned int> stab_dofs_rows2;
	stab_dofs_rows2.push_back(dof_indices_p[0]);
	system.matrix->add_matrix(Kstab2,stab_dofs_rows2,stab_dofs_cols2);
  test(4);
 
	//perf_log.pop("kstab");
  #endif 
  //endif USE_STAB