# Gwen's Handy Helpers
## A dirty list of instructions and links and references

### add a new package to a dev package's source
```
pgk() activate NaiveDynamics.jl
    add NewPackages

```
### getting a new REPL started after changing a type
I do believe using Revise is unnecessary, as VSCode supposedly already runs it. I havent noticed any issues, but I havent tested specifically if it works. Also, until I change the directory name, I must put NaiveDynamics.jl, not NaiveDynamics. That will mess life up!
```julia
pgk() dev ./NaiveDynamics.jl
using Revise
using NaiveDynamics
```
### TEST WITH LIVE LOCAL SERVER BEFORE DEPLOYING TO GITHUB AND HAVING TO CHAIN SMOKE COMMITS IN ORDER TO KNOW IF YOU
FORMATTED YOUR DOCS CORRECTLY