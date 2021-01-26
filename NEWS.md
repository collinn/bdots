# bdots 10-22-2020
- bdotsBoot successfully runs parallel on Windows
- TravisCI yaml file added
- test suite initialized
- added `bdRemove` for removing paired observations by fitCode
- no changes made to API
- data compressed in more efficient format
- package now satisfies *minimum* requirements for CRAN submission
- '...' no longer needed in custom curve functions

# bdots 10-31-2020
- Forgot to switch to devel branch, so pushing on master
- Now for real package passes requirements for CRAN submission
- generics and some documentation updated

# bdots 11-04-2020
- `p.adjust` exported to namespace to use `method == "oleson"`

# bdots 12-05-2020
- bdObj objects now have `fitCode` as integer
- polynomials currently implemented, though formating needs to be adjusted
- documentation updated for some functions
- plot functions in process of being disassembled, so expect them less cooperative than before
- Some steps towards curveList in bdotsBoot to having own class to be handled interactively more easily
