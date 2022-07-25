rm *.ll
for file in *.cpp;do clang++ -shared -fPIC $file -S -emit-llvm -o `echo $file| cut -d'.' -f1`.ll -O2 -fno-vectorize -fno-slp-vectorize -fno-unroll-loops -I/opt/python39/install/include/python3.9 -I/home/amit/Projects/COLABFIT/python_env/lib/python3.9/site-packages/pybind11/include -fPIC;done
for file in *.ll;do opt  $file  -load=/opt/enzyme/enzyme/build/Enzyme/LLVMEnzyme-12.so -enzyme -enzyme-cache-always=1 -o output_`echo $file| cut -d'.' -f1`.ll -S;done
clang++ -shared output_* -O3 -o descriptors`python3-config --extension-suffix` `python3-config --cflags --ldflags --libs `  -fPIC
 
