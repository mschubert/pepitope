.txdb_cache = new.env(parent=emptyenv())

.txdb_cache_reset = function(txdb) {
    rm(list=ls(envir=.txdb_cache, all.names=TRUE), envir=.txdb_cache)
    assign("txdb", txdb, envir=.txdb_cache)
    invisible(.txdb_cache)
}

.txdb_cache_current = function(txdb) {
    if (!exists("txdb", envir=.txdb_cache, inherits=FALSE) ||
            !identical(txdb, get("txdb", envir=.txdb_cache, inherits=FALSE))) {
        .txdb_cache_reset(txdb)
    }
    .txdb_cache
}

.txdb_cache_key = function(expr) {
    if (is.call(expr) && length(expr) > 1L)
        expr[[2L]] = quote(.txdb)
    deparse1(expr, width.cutoff=500L)
}

.txdb_cache_get = function(txdb, value) {
    key = .txdb_cache_key(substitute(value))
    cache = .txdb_cache_current(txdb)
    if (!exists(key, envir=cache, inherits=FALSE)) {
        assign(key, force(value), envir=cache)
    }
    get(key, envir=cache, inherits=FALSE)
}
