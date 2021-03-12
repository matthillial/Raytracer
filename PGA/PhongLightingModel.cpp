Color PhongLightingModel(){
		 float r = 0;
		 float g = 0;
		 float b = 0;
		 
		
		 //directional light calculations
			//Dir3D lightDir = dirdir * -1;
			//Dir3D normal = (hit - spherePos).normalized();
			//Dir3D viewDir = (eye - hit).normalized();
			//Dir3D halfway = (lightDir + viewDir).normalized();
		 
		 
		 //point light calculations
			Dir3D lightDir = (Point3D(plx,ply,plz) - hit).normalized();
	        Dir3D normal = (hit - spherePos).normalized();
	        Dir3D viewDir = (eye - hit).normalized();
			Dir3D halfway = (lightDir + viewDir).normalized();
			
		
		//cosine falloffs
		float difdot = (0 < dot(lightDir,normal)) ? dot(normal,lightDir) : 0;
		float specdot = (0 < pow(dot(halfway,normal),ns)) ? pow(dot(halfway,normal),ns) : 0;
		 
		 
		 //point light r g b values
			r = (ar * ambr) + (dr * (plr / pow(hit.distTo(Point3D(plx,ply,plz)),2)) * difdot) + (sr * (plr / pow(hit.distTo(Point3D(plx,ply,plz)),2)) * specdot);
			g = (ag * ambg) + (dg * (plg / pow(hit.distTo(Point3D(plx,ply,plz)),2)) * difdot) + (sg * (plg / pow(hit.distTo(Point3D(plx,ply,plz)),2)) * specdot);
			b = (ab * ambb) + (db * (plb / pow(hit.distTo(Point3D(plx,ply,plz)),2)) * difdot) + (sb * (plb / pow(hit.distTo(Point3D(plx,ply,plz)),2)) * specdot);
		 
		 //directional light r g b values
			//r = (ar * ambr) + (dr *  dirr  * difdot) + (sr * dirr * specdot);
			//g = (ag * ambg) + (dg *  dirg  * difdot) + (sg * dirg * specdot);
			//b = (ab * ambb) + (db *  dirb  * difdot) + (sb * dirb * specdot);
			
			
			
			return Color(r,g,b);
}
