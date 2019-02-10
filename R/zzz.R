# /////////////////////////////////////////////////////////////////////////
#
# Declare global variables and native symbol objects
#
# /////////////////////////////////////////////////////////////////////////

# Doing so, avoids the note from devtools::check():
# "no visible binding for global variable".
# See https://stackoverflow.com/a/12429344/5193830
# or https://stackoverflow.com/a/17807914/5193830

.onLoad <- function(...) {
  if (getRversion() >= "3.5.1")
    utils::globalVariables(
      c(
        'id', # solves - geobuffer_pts: no visible binding for global variable 'id'
        'as' # solves - .check_input: no visible global function definition for 'as'
      )
    )
}
