//The test probem is:
//Fix cube at the bottom in all directions.
//Then apply a force w.r.t the normal in the current configuration from the left.
//Fluid can flow out everywa

for (unsigned int s=0; s<elem->n_sides(); s++)
{
   if (elem->neighbor(s) == NULL)
   {   
			AutoPtr<Elem> side (elem->build_side(s));
			
			
	    Real hmax=(*elem).hmax();

			
			
  		fe_face_f->reinit(elem,s);  

      //fluid face 
  		const std::vector<std::vector<Real> >&  phi_face_f =  fe_face_f->get_phi();
  		const std::vector<std::vector<RealGradient> >& dphi_face_f = fe_face_f->get_dphi();
  		const std::vector<Real>& JxW_face_f = fe_face_f->get_JxW();
  		const std::vector<Point>& qface_point_f = fe_face_f->get_xyz();
  		const std::vector<Point>& face_normals_f = fe_face_f->get_normals();

			
			
			
			for (unsigned int qp=0; qp<qface_f->n_points(); qp++)
			{

					Real xf=qface_point_f[qp](0);
					Real yf=qface_point_f[qp](1);
					
					Point normal;
					normal(0)=face_normals_f[qp](0);
					normal(1)=face_normals_f[qp](1);

					#if THREED
					Real zf=qface_point_f[qp](2);
					normal(2)=face_normals_f[qp](2);
					#endif
										
					normal=normal.unit();
					
					Real pen_bc=DELTA_BC/hmax;
											
					for (unsigned int i=0; i<phi_face_f.size(); i++){			
	
						
						for (unsigned int j=0; j<phi_face_f.size(); j++){				
						
							
							#if !THREED
							//General
							Kuu(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_f[j][qp];  						
							
							Kvv(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_f[j][qp]; 
							#endif
							
							
							#if THREED
							//General
							Kuu(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(0)*phi_face_f[j][qp]*normal(0); 
							Kuv(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(1)*phi_face_f[j][qp]*normal(0);
							Kuw(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(2)*phi_face_f[j][qp]*normal(0); 
 						
							
							Kvu(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(0)*phi_face_f[j][qp]*normal(1); 
							Kvv(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(1)*phi_face_f[j][qp]*normal(1); 
							Kvw(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(2)*phi_face_f[j][qp]*normal(1); 
							
							
							Kwu(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(0)*phi_face_f[j][qp]*normal(2); 
							Kwv(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(1)*phi_face_f[j][qp]*normal(2); 
							Kvv(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(2)*phi_face_f[j][qp]*normal(2); 
							#endif
							
							
						}					
			}
			
			
		
					
			   		  
			
		} //end qp
		
		
	} //if (elem->neighbor(s) == NULL)
}// end boundary condition section  


