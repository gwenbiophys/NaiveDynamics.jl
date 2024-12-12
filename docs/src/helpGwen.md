# Gwen's Handy Helpers
### A dirty list of instructions and links and references

## add a new package to a dev package's source
```
pgk() activate NaiveDynamics.jl
    add NewPackages

```
## getting a new REPL started after changing a type
I do believe using Revise is unnecessary, as VSCode supposedly already runs it. I havent noticed any issues, but I havent tested specifically if it works. Also, until I change the directory name, I must put NaiveDynamics.jl, not NaiveDynamics. That will mess life up!
```julia
pgk() dev ./NaiveDynamics.jl
using Revise
using NaiveDynamics
```
## Deploy docs locally
- Get a local only make.jl with `deploydocs()` disabled. 
- Makesure git ignore includes the build directly entirely. 
- `serve(dir="/docs/build")` from the REPL. Have NaiveDynamics activated but make sure the local Julia environment has LiveServer. 
- If there is not yet a local Julia envionrment, type `activate` in pkg mode.

### After using 'markdown syntax' you have to press run on the markdown file
E.g. footnotes or links will not work until the .md itself is ran.
### ffmpeg for gif conversion, mayber tinker more for better qualty

`ffmpeg -t 4 -i /home/gwenk/Videos/iWant.mp4 -vf "fps=25,scale" -loop 1 /home/gwenk/Videos/iWant.gif`



## get a flame graph of allocations
by @profview_allocs
when just @profview will measure the execution time distribution

## a = b
if b is already defined, we are not asking julia to overwrite each value in a with the values in b, while maintaining a's structure, so long as a and b have the same structure.

We are asking Julia to ignore the value of a (if it has already been defined), and to point to b.

## To develop and use a package,
### either activate the REPL in the place you git cloned the repo
activate C:/Users/kucer/Desktop/julia/NaiveDynamics.jl
### or in the julia install itself where you may have to `dev` into the package.
activate C:/Users/Trist/.julia/dev/NaiveMD/NaiveDynamics.jl
using NaiveDynamics

## for type analysis
@code_warntype
