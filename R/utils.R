.pkg_state <- new.env(parent = emptyenv())
.pkg_state$nbda_msg_shown <- FALSE

.onAttach <- function(libname, pkgname) {
    if (!l10n_info()[["UTF-8"]]) {
        packageStartupMessage("Note: Non-UTF-8 locale detected. Some Unicode symbols may not display correctly.")
    }
}
