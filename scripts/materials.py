#! /usr/bin/python
#
# Materials: Writes the material cards for the Serpent input deck.


def write_materials(lib):
	'''Function to write material cards for Serpent input deck.
	Inputs: 
		lib:	String containing the neutron cross section library
				to use (e.g., 09c).
	Outputs:
		mats:	String containing the material cards'''
	
	mats = '''
%-------material definition--------------
%NOTE: VOLUMES OR MASS OF EACH MAY NEED TO 
%BE CALCULATED FOR BURNUP CALCULATIONS

%  FUEL SALT: 68.5-31.3-0.2 LiF-BeF2-UF4 at 900K 
%  DENSITY: 2.03434 g/cc
%  MELT TEMP: 450C or 742.15K
%  MATERIAL INFO FROM ONRL-4528 TALBE 3.1.
%  MAY NEED VOLUME OR MASS FOR BURNUP CALCULATIONS
mat fuel -2.03434  tmp 973
rgb 130 32 144
3006.{lib}   -0.000725     %  Li6
3007.{lib}  -14.495960     %  Li7
40090.{lib}  -8.508748     %  Be9
9019.{lib}  -75.588074     %  F19
92233.{lib}  -1.265844     %  U233
92234.{lib}  -0.140649     %  U234

%  BLANKET SALT: 71-27-2 LiF-ThF4-BeF2 at 900K
%  DENSITY: 4.43711 g/cc
%  MELT TEMP: 560C or 833.15K
%  MATERIAL INFO FROM ONRL-4528 TALBE 3.1.
%  MAY NEED VOLUME OR MASS FOR BURNUP CALCULATIONS
mat blanket -4.43711 tmp 973
rgb 0 157 254 
3006.{lib}   -0.000243     %  Li6
3007.{lib}   -4.855845     %  Li7
40090.{lib}  -0.175712     %  Be9
9019.{lib}  -33.892970     %  F19
90232.{lib} -61.052512     %  Th232
91233.{lib}  -0.022718     %  Pa233

%  NUCLEAR GRAPHITE: Natural concentration of carbon
%  DENSITY: 1.82 G/CC
mat graphite -1.82 moder graph 6000 tmp 973
rgb 130 130 130
6000.{lib} 1
%  THERMAL SCATTERING LIBRARY FOR GRAPHITE
therm graph gre7.08t

%  HELIUM: gas due to alpha particles
%  DENSITY: 54.19 E-6 g/cc
mat he -54.19E-6 tmp 973
rgb 0 255 16
2004.{lib} 1

% CONTROL ROD: NATURAL BORON at 900K
% DENSITY: 2.3 g/cc
% MELT TEMP: 2076C or 2349.15K
% 19.9 B10 and 80.1 B11
mat absorber -2.3 tmp 973
rgb 74 74 74
5010.{lib} -0.199
5011.{lib} -0.801

%  TODO: Hastelloy
mat hastelloy -8.86 tmp 973
rgb 139 69 19
28058.{lib}  -0.472120   %  Ni
 28060.{lib}  -0.181860   %  Ni
 28061.{lib}  -0.007905   %  Ni
 28062.{lib}  -0.025206   %  Ni
 28064.{lib}  -0.006419   %  Ni
 42100.{lib}  -0.015408   %  Mo
 42092.{lib}  -0.023744   %  Mo
 42094.{lib}  -0.014800   %  Mo
 42095.{lib}  -0.025472   %  Mo
 42096.{lib}  -0.026688   %  Mo
 42097.{lib}  -0.015280   %  Mo
 42098.{lib}  -0.038608   %  Mo
 24050.{lib}  -0.003041   %  Cr
 24052.{lib}  -0.058652   %  Cr
 24053.{lib}  -0.006651   %  Cr
 24054.{lib}  -0.001656   %  Cr
 26054.{lib}  -0.002923   %  Fe
 26056.{lib}  -0.045877   %  Fe
 26057.{lib}  -0.001059   %  Fe
 26058.{lib}  -0.000141   %  Fe
 14028.{lib}  -0.009223   %  Si
 14029.{lib}  -0.000468   %  Si
 14030.{lib}  -0.000309   %  Si
 25055.{lib}  -0.008000   %  Mn
 74182.{lib}  -0.001325   %  W
 74183.{lib}  -0.000715   %  W
 74184.{lib}  -0.001532   %  W
 74186.{lib}  -0.001422   %  W
 29063.{lib}  -0.002421   %  Cu
 29065.{lib}  -0.001079   %  Cu
'''
	mats = mats.format(**locals())

	return mats


if __name__ == '__main__':
	print "This is a module to write materials for the MSR core."
	raw_input("Press Ctrl+C to quit, or enter else to test it. ")
	print write_materials('09c')
