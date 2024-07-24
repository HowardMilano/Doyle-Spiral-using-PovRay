/**
 *
 * Based on Javascript code "Doyle spiral circle packing"
 * https://gist.github.com/robinhouston/6096562
 * Robin Houston
 *
 * Ported to PovRay and freely modified
 * July 23, 2024
 * Howard Milano
 *
 */

#include "functions.inc"
#include "rand.inc"
#include "shapes.inc"
#include "math.inc"

global_settings {
	assumed_gamma 2.2
	max_trace_level 100
}
 
camera { 
	location <0, 0, -90>
	direction <0.0, 0.0, 1.8>
	right <1.3, 0.0, 0.0>
	look_at <0.0, 0.0, 0.0>
}
 
#declare light_strength = 3.25; 
 
light_source {0*x color red light_strength green light_strength blue light_strength 
              translate <-60, 27, -80>}

#declare randomizer = seed(95235363);

#macro showOneSphere(object_radius, coloring)
object {
	#local  color_rgb = coloring / 256;
   	sphere {0, object_radius } 
  	texture {pigment {color rgbf <color_rgb, color_rgb, color_rgb>}
	    finish { phong 0.3 phong_size 60 } 
	}
}
#end

#macro dd(zz, tt, p, q)
	#local w = pow(zz, p/q);
	#local s = (p*tt + 2*pi) / q;
	(pow( zz*cos(tt) - w*cos(s), 2) + pow( zz*sin(tt) - w*sin(s), 2))
#end

#macro ddz_d(zz, tt, p, q) 
	#local w = pow(zz, p/q);
	#local s = (p*tt + 2*pi)/q;
	#local ddz_w = (p/q)*pow(zz, (p-q)/q);
	#local r1 = 	2*(w*cos(s) - zz*cos(tt))*(ddz_w*cos(s) - cos(tt));
	#local r2 = 2*(w*sin(s) - zz*sin(tt))*(ddz_w*sin(s) - sin(tt));
	(r1 + r2)
#end

#macro ddt_d(zz, tt, p, q)
	#local w = pow(zz, p/q);
	#local s = (p*tt + 2*pi)/q;
	#local dds_t = (p/q);
	#local r1 = 	2*( zz*cos(tt) - w*cos(s) )*( -zz*sin(tt) + w*sin(s)*dds_t );
	#local r2 = 2*( zz*sin(tt) - w*sin(s) )*( zz*cos(tt) - w*cos(s)*dds_t );
	(r1 + r2)
#end

#macro ss(zz, tt, p, q)
	pow(zz + pow(zz, p/q), 2)
#end

#macro ddz_s(zz, tt, p, q)
 	#local w = pow(zz, p/q);
	#local ddz_w = (p/q)*pow(zz, (p-q)/q);
	(2*(w+zz)*(ddz_w+1))
#end

#macro rr(zz, tt, p, q)
	(dd(zz,tt,p,q) / ss(zz,tt,p,q))
#end

#macro ddz_r(zz, tt, p, q)
	#local r1 = (ddz_d(zz,tt,p,q) * ss(zz,tt,p,q) - dd(zz,tt,p,q) * ddz_s(zz,tt,p,q));
	#local r2 = pow( ss(zz,tt,p,q), 2 ); 
	(r1 / r2)
#end

#macro ddt_r(zz, tt, p, q)
	((ddt_d(zz,tt,p,q) * ss(zz,tt,p,q)) / pow( ss(zz,tt,p,q), 2 ))
#end

#macro ff(zz, tt, p, q)
	(rr(zz,tt,0,1) - rr(zz,tt,p,q))
#end

#macro ddz_f(zz, tt, p, q)
	(ddz_r(zz,tt,0,1) - ddz_r(zz,tt,p,q))
#end

#macro ddt_f(zz, tt, p, q)
	(ddt_r(zz,tt,0,1) - ddt_r(zz,tt,p,q))
#end

#macro gg(zz, tt, p, q)
	#local gg1 = rr(zz,tt,0,1);
	#local gg2 = rr(pow(zz, p/q), (p*tt + 2*pi)/q, 0,1);
	(gg1 - gg2)
#end

#macro ddz_g(zz, tt, p, q)
	#local r1 = ddz_r(zz,tt,0,1);
	#local r2 = ddz_r(pow(zz, p/q), (p*tt + 2*pi)/q, 0,1);
	#local r3 = (p/q)*pow(zz, (p-q)/q);
	(r1 - r2 * r3)
#end

#macro ddt_g(zz, tt, p, q)
	(ddt_r(zz,tt,0,1) - ddt_r(pow(zz, p/q), (p*tt + 2*pi)/q, 0,1) * (p/q))
#end

#macro find_root(zz, tt, p, q, rootZ, rootT, rootR, rootSuccess)
	#local flag = 9;
	#local epsilon = 1e-10;
	
	// for(;;)
	#while (flag < 10) 
      	#local v_f = ff(zz, tt, p, q);
		#local v_g = gg(zz, tt, p, q);
		#if (-epsilon < v_f & v_f < epsilon & -epsilon < v_g & v_g < epsilon)
			#declare rootZ = zz;
			#declare rootT = tt;
			#declare rootR = sqrt(rr(zz,tt,0,1));
			#declare rootSuccess = 1;
			#local flag = 10;
			// #debug "Found root\n"
		#else
			#local a = ddz_f(zz,tt, p, q);
                	#local b = ddt_f(zz,tt, p, q);
                	#local c = ddz_g(zz,tt, p, q);
                	#local d = ddt_g(zz,tt, p, q);
                	#local det = a*d-b*c;
			#if (-epsilon < det & det < epsilon)
				#local flag = 10;
				#declare rootSuccess = 0;
				#debug concat("No root for det ", str(det, 0, 5), "\n")
			#else
				#local zz = zz - (d*v_f - b*v_g)/det;
				#local tt = tt - (a*v_g - c*v_f)/det;
				#if (zz < epsilon)
					#local flag = 10;
					#declare rootSuccess = 0;
					#debug concat("No root for zz ", str(zz, 0, 5), "\n")
				#end
			#end
		#end
	#end
#end

#macro doyle(p, q, doyleA, doyleB, doyleR, doyleModa, doyleArga,  doyleSuccess)
	#declare rootZ = -1;
	#declare rootT = -1;
	#declare rootR = -1;
	#declare rootSuccess = -1;
	find_root(2, 0, p, q, rootZ, rootT, rootR, rootSuccess)
      #if (rootSuccess < 1)
      	#debug "Failed to find root\n"
      	#declare doyleSuccess = 0;
      #else
		#declare doyleA = <rootZ * cos(rootT), rootZ * sin(rootT)>;
		#local corootZ = pow(rootZ, p/q);
		#local corootT = (p*rootT+2*pi)/q;
		#declare doyleB = <corootZ * cos(corootT), corootZ * sin(corootT)>;
		#declare doyleR = rootR;
		#declare doyleModa = rootZ;
		#declare doyleArga = rootT;
		#declare doyleSuccess = 1;
	#end
#end

#macro cmul (ww, zz, cmulOut)
	#declare cmulOut =  <ww.u*zz.u - ww.v*zz.v, ww.u*zz.v + ww.v*zz.u>;
#end

#macro modulus(p)
	sqrt(p.u*p.u + p.v*p.v)
#end

#macro crecip(delta, crecipOut)
	#local d = delta.u*delta.u + delta.v*delta.v;
	#declare crecipOut = <delta.u/d, -delta.v/d>;
#end

#macro doyle_spiral(spiral_radius, startPoint, delta, minD, maxD)
	#local ball_color = rand(randomizer) * 64 + 64;
	#declare crecipOut = <0, 0>;
	crecip(delta, crecipOut)
	
	#local mod_delta = modulus(delta);
	#local mod_recip_delta = 1/mod_delta;
        
      // Spiral outwards
      #local qq = startPoint;
      #local mod_q = modulus(qq);
      #while (mod_q < maxD)
      	// #local str1 = str(qq.u, 0, 5);
      	// #local str2 = str(qq.v, 0, 5);
      	// #local str3 = str(mod_q*spiral_radius, 0, 5);
      	// #debug concat("circle at ", str1, ", ", str2, ", r = ", str3, "\n")
      	object {
	      	showOneSphere(mod_q * spiral_radius, ball_color)
	      	translate<qq.u, qq.v, 0>
      	}
      	#declare cmulOut = <0, 0>;
      	cmul(qq, delta, cmulOut)
      	#local qq = cmulOut;
      	#local mod_q = mod_q * mod_delta;
      #end

	// Spiral inwards
      #declare cmulOut = <0, 0>;
      cmul(startPoint, crecipOut, cmulOut)
      #local qq = cmulOut;
      #local mod_q = modulus(qq);
      #while (mod_q > min_d)
      	// #local str1 = str(qq.u, 0, 5);
      	// #local str2 = str(qq.v, 0, 5);
      	// #local str3 = str(mod_q*spiral_radius, 0, 5);
      	// #debug concat("circle at ", str1, ", ", str2, ", r = ", str3, "\n")
      	object {
	      	showOneSphere(mod_q * spiral_radius, ball_color)
	      	translate<qq.u, qq.v, 0>
      	}
		#declare cmulOut = <0, 0>;
		cmul(qq, crecipOut, cmulOut)
		#local qq = cmulOut;
		#local mod_q = mod_q * mod_recip_delta;        	
	#end

#end

//
// Main Program
//

#local p = 4;
#local q = 29;

#declare doyleA = <0, 0>;
#declare doyleB = <0, 0>;
#declare doyleR = 0;
#declare doyleModa = 0;
#declare doyleArga = 0;
#declare doyleSuccess = 0;	
doyle(p, q, doyleA, doyleB, doyleR, doyleModa, doyleArga,  doyleSuccess)
#if (doyleSuccess < 1)
	#debug "Calling doyle failed\n"
#else
	union {
		#local min_d = 0.1;
		#local max_d = 45;
		#local start = doyleA;
		#local spiralCount = 0;
		#while (spiralCount < q)
			doyle_spiral(doyleR, start, doyleA, min_d, max_d)
			#declare cmulOut = <0, 0>;
			cmul(start, doyleB, cmulOut)
			#local start = cmulOut;
			#local spiralCount = spiralCount + 1;
		#end
		scale<1, 1, 1>
		translate<0, 0, 0>
	}
#end
