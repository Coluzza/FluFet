CLAPACK=$(CLAPACK_LIB)
LIB=lib
INCLUDE=-I$(CLAPACK_INCLUDE)
CFLAGS=  -O3 -g3 -fno-omit-frame-pointer -gdwarf-4


all:

	gcc $(CFLAGS) -c -static contactmap-Shae_NY_design-2.c -o contactmap-Shae_NY_design-2.o
	gcc $(CFLAGS) -c -static contactmap-Shae_NY_design-GO.c -o contactmap-Shae_NY_design-GO.o
	gcc $(CFLAGS) -c -static contactmap-Shae_NY_design-2.c -o contactmap-Shae_NY_design-2_GO.o -DGO
	gcc $(CFLAGS) -c -static contactmap-Shae_NY_design-Graphene.c -o contactmap-Shae_NY_design-Graphene.o -DGO
	gcc $(CFLAGS) -c -static contactmap-Shae_NY_design-COVID.c -o contactmap-Shae_NY_design-COVID.o -DGO
	gcc $(CFLAGS) -c -static contactmap-Shae_NY_design-NanoPart.c -o contactmap-Shae_NY_design-NanoPart.o -DGO
	gcc $(CFLAGS) -c -static contactmap-Shae_NY_design-NanoPart_multy.c -o contactmap-Shae_NY_design-NanoPart_multy.o -DGO
	gcc $(CFLAGS) -c -static contactmap-Shae_NY_design-GraftFlow.c -o contactmap-Shae_NY_design-GraftFlow.o -DGO
	gcc $(CFLAGS) -c -static contactmap-Shae_NY_design-VirusFlow.c -o contactmap-Shae_NY_design-VirusFlow.o -DGO
	gcc $(CFLAGS) -c -static contactmap-Shae_NY_design-2.c -o contactmap-Shae_NY_design-2_GAUSS.o -DGAUSS -DGO
	gcc $(CFLAGS) -c -static GraftFlow_Brush.c -o GraftFlow_Brush.o
	
	gcc $(CFLAGS) -c -static Quantasome.c -o Quantasome.o
	gcc $(CFLAGS) -c -static Quantasome_charges.c -o Quantasome_charges.o
	gcc $(CFLAGS) -c $(LIB)/sort2.c -o $(LIB)/sort2.o
	gcc $(CFLAGS) -c $(LIB)/partint.c -o $(LIB)/partint.o
	gcc $(CFLAGS) -c $(LIB)/ran3.c -o $(LIB)/ran3.o
	gcc $(CFLAGS) -c $(LIB)/param_parser.c -o $(LIB)/param_parser.o
	
	gcc -Wall -Wno-unused-variable $(CFLAGS) contactmap-Shae_NY_design-2.o $(LIB)/sort2.o $(LIB)/partint.o  -o contactmap-Shae_NY_design-2.x  -lm   	
	gcc -Wall -Wno-unused-variable $(CFLAGS) contactmap-Shae_NY_design-2_GO.o $(LIB)/sort2.o $(LIB)/partint.o  -o contactmap-Shae_NY_design-2_GO.x  -lm  
	gcc -Wall -Wno-unused-variable $(CFLAGS) contactmap-Shae_NY_design-GO.o $(LIB)/sort2.o $(LIB)/partint.o $(LIB)/param_parser.o -o contactmap-Shae_NY_design-GO.x  -lm  
	gcc -Wall -Wno-unused-variable $(CFLAGS) contactmap-Shae_NY_design-Graphene.o $(LIB)/sort2.o $(LIB)/partint.o  -o contactmap-Shae_NY_design-Graphene.x  -lm  
	gcc -Wall -Wno-unused-variable $(CFLAGS) contactmap-Shae_NY_design-COVID.o $(LIB)/sort2.o $(LIB)/partint.o  $(LIB)/param_parser.o  -o contactmap-Shae_NY_design-COVID.x  -lm   
	gcc -Wall -Wno-unused-variable $(CFLAGS) contactmap-Shae_NY_design-2_GAUSS.o $(LIB)/sort2.o $(LIB)/partint.o  -o contactmap-Shae_NY_design-2_GAUSS.x  -lm   
	gcc -Wall -Wno-unused-variable $(CFLAGS) contactmap-Shae_NY_design-NanoPart.o $(LIB)/sort2.o $(LIB)/partint.o $(LIB)/ran3.o $(LIB)/param_parser.o -o contactmap-Shae_NY_design-NanoPart.x  -lm
	gcc -Wall -Wno-unused-variable $(CFLAGS) contactmap-Shae_NY_design-NanoPart_multy.o $(LIB)/sort2.o $(LIB)/partint.o $(LIB)/ran3.o $(LIB)/param_parser.o -o contactmap-Shae_NY_design-NanoPart_multy.x  -lm
	gcc -Wall -Wno-unused-variable $(CFLAGS) contactmap-Shae_NY_design-GraftFlow.o $(LIB)/sort2.o $(LIB)/partint.o $(LIB)/ran3.o $(LIB)/param_parser.o -o contactmap-Shae_NY_design-GraftFlow.x  -lm
	gcc -Wall -Wno-unused-variable $(CFLAGS) contactmap-Shae_NY_design-VirusFlow.o $(LIB)/sort2.o $(LIB)/partint.o $(LIB)/ran3.o $(LIB)/param_parser.o -o contactmap-Shae_NY_design-VirusFlow.x  -lm
	gcc -Wall -Wno-unused-variable $(CFLAGS) GraftFlow_Brush.o  $(LIB)/sort2.o $(LIB)/partint.o $(LIB)/ran3.o $(LIB)/param_parser.o  -o GraftFlow_Brush.x  -lm	
	gcc -Wall -Wno-unused-variable $(CFLAGS) Quantasome.o $(LIB)/param_parser.o -o Quantasome.x  -lm
	gcc -Wall -Wno-unused-variable $(CFLAGS) Quantasome_charges.o $(LIB)/param_parser.o -o Quantasome_charges.x  -lm


general:
	gcc $(CFLAGS) -c -static GraftFlow_Brush_general.c -o GraftFlow_Brush_general.o
	gcc $(CFLAGS) -c $(LIB)/sort2.c -o $(LIB)/sort2.o
	gcc $(CFLAGS) -c $(LIB)/partint.c -o $(LIB)/partint.o
	gcc $(CFLAGS) -c $(LIB)/ran3.c -o $(LIB)/ran3.o
	gcc $(CFLAGS) -c $(LIB)/param_parser.c -o $(LIB)/param_parser.o
	gcc -Wall -Wno-unused-variable $(CFLAGS) GraftFlow_Brush_general.o  $(LIB)/sort2.o $(LIB)/partint.o $(LIB)/ran3.o $(LIB)/param_parser.o  -o GraftFlow_Brush_general.x  -lm	
