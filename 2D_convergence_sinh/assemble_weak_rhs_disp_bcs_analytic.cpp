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

				
				Number value_u = exact_2D_solution_u(qface_point_f[qp], es.parameters,"null","void");
				Number value_v = exact_2D_solution_v(qface_point_f[qp], es.parameters,"null","void");
            
				Real xf=qface_point_f[qp](0);
				Real yf=qface_point_f[qp](1);

						
				Point disp_val;
				disp_val(0)=value_u;
				disp_val(1)=value_v;

						
				Point normal;
				normal(0)=face_normals_f[qp](0);
				normal(1)=face_normals_f[qp](1);
				
				
				#if THREED
				normal(2)=face_normals_f[qp](2);
				Real zf=qface_point_f[qp](2);
				Number value_z = exact_2D_solution_z(qface_point_f[qp], es.parameters,"null","void");
				disp_val(2)=value_z;
				#endif
				
				normal=normal.unit();
	
	
				Real pen_bc=DELTA_BC/hmax;
					
				for (unsigned int i=0; i<phi_face_f.size(); i++){			
		 		
#if !THREED
					Fu(i) +=  pen_bc*JxW_face_f[qp]*disp_val(0)*phi_face_f[i][qp];
 						
				  Fv(i) +=  pen_bc*JxW_face_f[qp]*disp_val(1)*phi_face_f[i][qp];
#endif
					
					
#if THREED
					Fu(i) +=  pen_bc*JxW_face_f[qp]*disp_val(0)*normal(0)*phi_face_f[i][qp]*normal(0);
				  Fu(i) +=  pen_bc*JxW_face_f[qp]*disp_val(1)*normal(1)*phi_face_f[i][qp]*normal(0);
					Fu(i) +=  pen_bc*JxW_face_f[qp]*disp_val(2)*normal(2)*phi_face_f[i][qp]*normal(0);
 						
				  Fv(i) +=  pen_bc*JxW_face_f[qp]*disp_val(0)*normal(0)*phi_face_f[i][qp]*normal(1);
				  Fv(i) +=  pen_bc*JxW_face_f[qp]*disp_val(1)*normal(1)*phi_face_f[i][qp]*normal(1);
					Fv(i) +=  pen_bc*JxW_face_f[qp]*disp_val(2)*normal(2)*phi_face_f[i][qp]*normal(1);
					
					
					Fw(i) +=  pen_bc*JxW_face_f[qp]*disp_val(0)*normal(0)*phi_face_f[i][qp]*normal(2);
				  Fw(i) +=  pen_bc*JxW_face_f[qp]*disp_val(1)*normal(1)*phi_face_f[i][qp]*normal(2);
					Fw(i) +=  pen_bc*JxW_face_f[qp]*disp_val(2)*normal(2)*phi_face_f[i][qp]*normal(2);
#endif
				}
				

		} //end qp
	} //if (elem->neighbor(s) == NULL)
}// end boundary condition section  


