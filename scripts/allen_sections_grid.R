# NOTE:
#   The Allen adult mouse expression grid is 200 µm (67 x 41 x 58).
#   This backend is appropriate for voxelized signatures and coarse single-gene plots.
#   For sharp single-gene ISH panels, use SectionImage instead of grid_data.
#
# Tested logic: synthetic geometry and internal function flow.
# Live Allen downloads require internet access in the user's R session.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Please install 'jsonlite': install.packages('jsonlite')")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

.ALLEN_CACHE_DIR <- file.path(tempdir(), "allen_grid_cache_v11")
if (!dir.exists(.ALLEN_CACHE_DIR)) dir.create(.ALLEN_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

clear_allen_cache <- function() {
  if (dir.exists(.ALLEN_CACHE_DIR)) unlink(.ALLEN_CACHE_DIR, recursive = TRUE, force = TRUE)
  dir.create(.ALLEN_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  invisible(.ALLEN_CACHE_DIR)
}

# ----------------------------
# URL / download helpers
# ----------------------------

build_url <- function(base_url, query = list()) {
  if (length(query) == 0) return(base_url)
  q <- vapply(names(query), function(k) {
    val <- query[[k]]
    if (length(val) > 1) val <- paste(val, collapse = ",")
    paste0(utils::URLencode(k, reserved = TRUE), "=",
           utils::URLencode(as.character(val), reserved = TRUE))
  }, character(1))
  paste0(base_url, if (grepl("\\?", base_url, fixed = TRUE)) "&" else "?", paste(q, collapse = "&"))
}

download_text <- function(url) {
  paste(readLines(url, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

download_file_retry <- function(url, dest, quiet = TRUE) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest)) unlink(dest, force = TRUE)
  utils::download.file(url, destfile = dest, mode = "wb", quiet = quiet)
  if (!file.exists(dest) || is.na(file.info(dest)$size) || file.info(dest)$size <= 0) {
    stop("Download failed or returned empty file: ", url)
  }
  dest
}

api_get_json <- function(url, query = list(), simplify = FALSE) {
  full_url <- build_url(url, query)
  txt <- tryCatch(download_text(full_url), error = function(e) {
    stop("Failed to retrieve:\n", full_url, "\n", conditionMessage(e))
  })
  out <- jsonlite::fromJSON(txt, simplifyDataFrame = simplify)
  success_value <- out$success %||% TRUE
  success_flag <- isTRUE(success_value) || identical(tolower(as.character(success_value)[1]), "true")
  if (!success_flag) {
    err <- out$msg %||% out$message %||% txt
    stop("Allen API returned success = FALSE. Details: ", paste(err, collapse = " "))
  }
  out$msg
}

read_remote_csv <- function(url) {
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)
  download_file_retry(url, tf)
  utils::read.csv(tf, stringsAsFactors = FALSE, check.names = FALSE)
}

# ----------------------------
# Allen dataset lookup
# ----------------------------

normalize_records <- function(x) {
  if (is.null(x)) return(list())
  if (is.data.frame(x)) return(split(x, seq_len(nrow(x))))
  if (!is.list(x)) return(list())
  if (length(x) == 0) return(list())
  if (!is.null(names(x)) && any(nzchar(names(x)))) return(list(x))
  x
}

pluck_chr <- function(x, path, default = "") {
  cur <- x
  for (nm in path) {
    if (is.null(cur) || is.null(cur[[nm]])) return(default)
    cur <- cur[[nm]]
  }
  if (is.null(cur) || length(cur) == 0) return(default)
  as.character(cur[[1]])
}

collapse_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return("")
  paste(unique(x), collapse = ";")
}

as_flag <- function(x) {
  if (is.null(x) || length(x) == 0) return(FALSE)
  if (is.logical(x)) return(isTRUE(x[1]))
  tolower(as.character(x[1])) %in% c("true", "t", "1", "yes", "y")
}

parse_dataset_records <- function(raw_msg, fallback_gene = NA_character_) {
  if (is.null(raw_msg) || length(raw_msg) == 0) {
    return(data.frame(
      id = integer(), gene = character(), plane = character(),
      probe_orientation = character(), treatments = character(),
      delegate = logical(), stringsAsFactors = FALSE
    ))
  }

  recs <- if (is.data.frame(raw_msg)) {
    split(raw_msg, seq_len(nrow(raw_msg)))
  } else if (is.list(raw_msg) && length(raw_msg) > 0 &&
             all(vapply(raw_msg, is.list, logical(1)))) {
    raw_msg
  } else {
    list(raw_msg)
  }

  rows <- lapply(recs, function(x) {
    probes <- normalize_records(x$probes)
    treatments <- normalize_records(x$treatments)
    genes <- normalize_records(x$genes)

    probe_orientations <- vapply(probes, function(p) {
      # orientation may be direct string or nested list/data.frame
      if (is.list(p$orientation) || is.data.frame(p$orientation)) {
        tolower(pluck_chr(p, c("orientation", "name"), default = ""))
      } else {
        val <- p$orientation %||% ""
        tolower(as.character(val)[1])
      }
    }, character(1))

    treatment_names <- vapply(treatments, function(t) {
      tolower(pluck_chr(t, c("name"), default = ""))
    }, character(1))

    gene_acronyms <- vapply(genes, function(g) {
      pluck_chr(g, c("acronym"), default = "")
    }, character(1))

    data.frame(
      id = suppressWarnings(as.integer(x$id %||% NA_integer_)),
      gene = if (length(gene_acronyms) > 0 && nzchar(gene_acronyms[1])) gene_acronyms[1] else fallback_gene,
      plane = tolower(pluck_chr(x, c("plane_of_section", "name"), default = "")),
      probe_orientation = collapse_chr(probe_orientations),
      treatments = collapse_chr(treatment_names),
      delegate = as_flag(x$delegate),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[!is.na(out$id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

parse_dataset_csv <- function(df, fallback_gene = NA_character_) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      id = integer(), gene = character(), plane = character(),
      probe_orientation = character(), treatments = character(),
      delegate = logical(), stringsAsFactors = FALSE
    ))
  }

  nm <- tolower(names(df))
  names(df) <- nm

  pick <- function(cands) {
    hit <- intersect(tolower(cands), nm)
    if (length(hit) == 0) return(rep(NA_character_, nrow(df)))
    as.character(df[[hit[1]]])
  }

  out <- data.frame(
    id = suppressWarnings(as.integer(pick(c("id", "section_data_set_id", "data_sets.id")))),
    gene = pick(c("gene", "genes.acronym")),
    plane = tolower(pick(c("plane", "plane_of_sections.name"))),
    probe_orientation = tolower(pick(c("probe_orientation", "probes.orientation"))),
    treatments = tolower(pick(c("treatments", "treatment", "treatments.name"))),
    delegate = tolower(pick(c("delegate", "data_sets.delegate"))) %in% c("true", "t", "1", "yes"),
    stringsAsFactors = FALSE
  )
  bad <- is.na(out$gene) | !nzchar(out$gene)
  out$gene[bad] <- fallback_gene
  out <- out[!is.na(out$id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

fetch_section_datasets_raw <- function(gene, product = "Mouse", verbose = FALSE) {
  gene <- as.character(gene)[1]

  crit_direct <- paste0(
    "[failed$eq'false'],",
    "products[abbreviation$eq'", product, "'],",
    "genes[acronym$eq'", gene, "']"
  )

  crit_rma <- paste0(
    "model::SectionDataSet,rma::criteria,",
    "[failed$eq'false'],",
    "products[abbreviation$eq'", product, "'],",
    "genes[acronym$eq'", gene, "']"
  )

  crit_csv <- paste0(
    "model::SectionDataSet,rma::criteria,",
    "[failed$eq'false'],",
    "products[abbreviation$eq'", product, "'],",
    "genes[acronym$eq'", gene, "'],",
    "genes,plane_of_section,treatments,probes,",
    "rma::options,",
    "[tabular$eq'data_sets.id+as+id','genes.acronym+as+gene',",
    "'plane_of_sections.name+as+plane','probes.orientation+as+probe_orientation',",
    "'treatments.name+as+treatments','data_sets.delegate+as+delegate'],",
    "[order$eq'data_sets.id']"
  )

  attempts <- list(
    list(
      name = "direct_json",
      fun = function() api_get_json(
        "https://api.brain-map.org/api/v2/data/SectionDataSet/query.json",
        query = list(criteria = crit_direct,
                     include = "genes,plane_of_section,probes,treatments",
                     num_rows = "all"),
        simplify = FALSE
      ),
      parser = function(x) parse_dataset_records(x, fallback_gene = gene)
    ),
    list(
      name = "rma_json",
      fun = function() api_get_json(
        "https://api.brain-map.org/api/v2/data/query.json",
        query = list(criteria = crit_rma,
                     include = "genes,plane_of_section,probes,treatments",
                     num_rows = "all"),
        simplify = FALSE
      ),
      parser = function(x) parse_dataset_records(x, fallback_gene = gene)
    ),
    list(
      name = "csv_tabular",
      fun = function() {
        url <- paste0(
          "https://api.brain-map.org/api/v2/data/query.csv?criteria=",
          utils::URLencode(crit_csv, reserved = TRUE),
          "&num_rows=all"
        )
        read_remote_csv(url)
      },
      parser = function(x) parse_dataset_csv(x, fallback_gene = gene)
    )
  )

  for (att in attempts) {
    res <- tryCatch(att$fun(), error = function(e) e)
    if (inherits(res, "error")) {
      if (verbose) message("Lookup attempt ", att$name, " failed for ", gene, ": ", conditionMessage(res))
      next
    }
    parsed <- tryCatch(att$parser(res), error = function(e) e)
    if (inherits(parsed, "error")) {
      if (verbose) message("Parse attempt ", att$name, " failed for ", gene, ": ", conditionMessage(parsed))
      next
    }
    if (nrow(parsed) > 0) {
      rownames(parsed) <- NULL
      return(parsed)
    }
    if (verbose) message("Lookup attempt ", att$name, " returned 0 rows for ", gene)
  }

  data.frame(
    id = integer(), gene = character(), plane = character(),
    probe_orientation = character(), treatments = character(),
    delegate = logical(), stringsAsFactors = FALSE
  )
}

inspect_allen_datasets <- function(gene,
                                   plane = NULL,
                                   orientation = NULL,
                                   product = "Mouse",
                                   require_ish = FALSE,
                                   verbose = FALSE) {
  ds <- fetch_section_datasets_raw(gene = gene, product = product, verbose = verbose)

  if (!is.null(plane) && nrow(ds) > 0 && any(nzchar(ds$plane))) {
    ds <- ds[tolower(ds$plane) == tolower(plane), , drop = FALSE]
  }

  if (!is.null(orientation) && nrow(ds) > 0 && any(nzchar(ds$probe_orientation))) {
    keep <- !nzchar(ds$probe_orientation) | grepl(tolower(orientation), ds$probe_orientation, fixed = TRUE)
    tmp <- ds[keep, , drop = FALSE]
    if (nrow(tmp) > 0) ds <- tmp
  }

  if (require_ish && nrow(ds) > 0 && any(nzchar(ds$treatments))) {
    tmp <- ds[grepl("ish", ds$treatments, ignore.case = TRUE), , drop = FALSE]
    if (nrow(tmp) > 0) ds <- tmp
  }

  if (nrow(ds) == 0) return(ds)
  ds <- ds[!duplicated(ds$id), , drop = FALSE]
  ds <- ds[order(ds$plane, ds$id), , drop = FALSE]
  rownames(ds) <- NULL
  ds
}

find_section_datasets <- function(gene,
                                  plane = c("coronal", "sagittal"),
                                  orientation = "antisense",
                                  product = "Mouse",
                                  require_ish = TRUE,
                                  verbose = FALSE) {
  plane <- match.arg(plane)
  ds <- inspect_allen_datasets(
    gene = gene,
    plane = plane,
    orientation = orientation,
    product = product,
    require_ish = require_ish,
    verbose = verbose
  )

  if (nrow(ds) == 0) {
    stop("No experiments found for gene ", gene, " (plane = ", plane, ")")
  }
  ds
}

# ----------------------------
# MHD / RAW reading
# ----------------------------

parse_mhd <- function(mhd_path) {
  lines <- readLines(mhd_path, warn = FALSE)
  kv <- strsplit(lines, "=", fixed = TRUE)
  kv <- kv[lengths(kv) == 2]
  keys <- trimws(vapply(kv, `[`, character(1), 1))
  vals <- trimws(vapply(kv, `[`, character(1), 2))
  stats::setNames(as.list(vals), keys)
}

read_uint32_le <- function(con, n) {
  b <- readBin(con, "integer", n = 4L * n, size = 1L, signed = FALSE, endian = "little")
  m <- matrix(b, nrow = 4L)
  out <- m[1, ] + 256 * m[2, ] + 65536 * m[3, ] + 16777216 * m[4, ]
  as.numeric(out)
}

read_mhd <- function(mhd_path) {
  hdr <- parse_mhd(mhd_path)

  dims <- as.integer(strsplit(hdr[["DimSize"]], "\\s+")[[1]])
  etype <- hdr[["ElementType"]]
  data_file <- hdr[["ElementDataFile"]]

  msb <- hdr[["ElementByteOrderMSB"]]
  if (is.null(msb) || is.na(msb) || msb == "") {
    msb <- hdr[["BinaryDataByteOrderMSB"]]
  }
  endian <- if (is.null(msb) || is.na(msb) || msb == "" || tolower(msb) %in% c("false", "0")) "little" else "big"

  raw_path <- file.path(dirname(mhd_path), data_file)
  con <- file(raw_path, "rb")
  on.exit(close(con), add = TRUE)

  n <- prod(dims)
  vec <- switch(
    etype,
    MET_FLOAT = readBin(con, what = numeric(), n = n, size = 4L, endian = endian),
    MET_DOUBLE = readBin(con, what = numeric(), n = n, size = 8L, endian = endian),
    MET_SHORT = readBin(con, what = integer(), n = n, size = 2L, signed = TRUE, endian = endian),
    MET_USHORT = readBin(con, what = integer(), n = n, size = 2L, signed = FALSE, endian = endian),
    MET_UCHAR = readBin(con, what = integer(), n = n, size = 1L, signed = FALSE, endian = endian),
    MET_CHAR = readBin(con, what = integer(), n = n, size = 1L, signed = TRUE, endian = endian),
    MET_INT = readBin(con, what = integer(), n = n, size = 4L, signed = TRUE, endian = endian),
    MET_UINT = {
      if (endian != "little") stop("MET_UINT currently implemented only for little-endian data")
      read_uint32_le(con, n)
    },
    stop("Unsupported ElementType: ", etype)
  )

  if (length(vec) != n) {
    stop("Failed to read full volume from ", raw_path)
  }
  array(vec, dim = dims)
}

# ----------------------------
# Downloads: expression + annotation
# ----------------------------

download_expression_grid <- function(dataset_id, metric = c("energy", "density", "intensity"), verbose = FALSE) {
  metric <- match.arg(metric)
  outdir <- file.path(.ALLEN_CACHE_DIR, paste0("grid_", dataset_id, "_", metric))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  zipfile <- file.path(outdir, paste0(dataset_id, "_", metric, ".zip"))
  mhd <- list.files(outdir, pattern = paste0("^", metric, ".*\\.mhd$|", metric, "\\.mhd$"),
                    full.names = TRUE, recursive = TRUE)

  if (length(mhd) == 0) {
    url <- paste0("https://api.brain-map.org/grid_data/download/", dataset_id, "?include=", metric)
    if (verbose) message("Downloading grid: ", url)
    download_file_retry(url, zipfile)
    utils::unzip(zipfile, exdir = outdir)
    mhd <- list.files(outdir, pattern = paste0("^", metric, ".*\\.mhd$|", metric, "\\.mhd$"),
                      full.names = TRUE, recursive = TRUE)
  }

  if (length(mhd) == 0) {
    stop("Could not find MHD file for metric ", metric, " in dataset ", dataset_id)
  }

  vol <- read_mhd(mhd[1])
  vol[is.finite(vol) & vol < 0] <- NA_real_
  vol
}

download_annotation_grid <- function(verbose = FALSE) {
  outdir <- file.path(.ALLEN_CACHE_DIR, "annotation_200um")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  mhd <- list.files(outdir, pattern = "gridAnnotation.*\\.mhd$|annotation.*\\.mhd$|\\.mhd$",
                    full.names = TRUE, recursive = TRUE)

  if (length(mhd) == 0) {
    urls <- c(
      "https://download.alleninstitute.org/informatics-archive/current-release/mouse_annotation/P56_Mouse_gridAnnotation.zip",
      "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/P56_Mouse_gridAnnotation.zip"
    )
    zipfile <- file.path(outdir, "P56_Mouse_gridAnnotation.zip")
    ok <- FALSE
    for (u in urls) {
      ok <- tryCatch({
        if (verbose) message("Trying annotation URL: ", u)
        download_file_retry(u, zipfile)
        TRUE
      }, error = function(e) FALSE)
      if (ok) break
    }
    if (!ok) stop("Failed to download P56_Mouse_gridAnnotation.zip from official Allen URLs")
    utils::unzip(zipfile, exdir = outdir)
    mhd <- list.files(outdir, pattern = "gridAnnotation.*\\.mhd$|annotation.*\\.mhd$|\\.mhd$",
                      full.names = TRUE, recursive = TRUE)
  }

  if (length(mhd) == 0) stop("Could not find annotation MHD file after unzip")

  ann <- read_mhd(mhd[1])
  ann
}

# ----------------------------
# Structure graph
# ----------------------------

flatten_structure_tree <- function(node, parent_id = NA_integer_, depth = 0L, out = NULL) {
  row <- data.frame(
    id = as.integer(node$id),
    parent_structure_id = as.integer(parent_id),
    depth = as.integer(depth),
    acronym = as.character(node$acronym %||% ""),
    name = as.character(node$name %||% ""),
    color = paste0("#", as.character(node$color_hex_triplet %||% "DDDDDD")),
    stringsAsFactors = FALSE
  )
  out <- if (is.null(out)) row else rbind(out, row)
  children <- node$children %||% list()
  if (length(children) > 0) {
    for (ch in children) {
      out <- flatten_structure_tree(ch, parent_id = as.integer(node$id), depth = depth + 1L, out = out)
    }
  }
  out
}

download_structure_graph <- function(graph_id = 1L, verbose = FALSE) {
  cache_file <- file.path(.ALLEN_CACHE_DIR, paste0("structure_graph_", graph_id, ".json"))
  if (!file.exists(cache_file)) {
    url <- paste0("https://api.brain-map.org/api/v2/structure_graph_download/", graph_id, ".json")
    if (verbose) message("Downloading structure graph: ", url)
    txt <- download_text(url)
    writeLines(txt, cache_file)
  }
  js <- jsonlite::fromJSON(cache_file, simplifyDataFrame = FALSE)
  msg <- js$msg
  if (is.null(msg) || length(msg) == 0) stop("Structure graph JSON did not contain msg")
  tree <- msg[[1]]
  out <- flatten_structure_tree(tree)
  rownames(out) <- NULL
  out
}

collapse_ids_to_depth <- function(id_mat, structure_df, atlas_depth = 5L) {
  ids <- structure_df$id
  parent <- structure_df$parent_structure_id
  depth <- structure_df$depth
  names(parent) <- as.character(ids)
  names(depth) <- as.character(ids)

  collapse_one <- function(id) {
    if (is.na(id) || id <= 0) return(NA_integer_)
    cur <- as.character(as.integer(id))
    guard <- 0L
    while (guard < 100L) {
      d <- depth[[cur]]
      p <- parent[[cur]]
      if (is.null(d) || is.na(d)) return(as.integer(id))
      if (d <= atlas_depth || is.null(p) || is.na(p)) return(as.integer(cur))
      cur <- as.character(as.integer(p))
      guard <- guard + 1L
    }
    as.integer(id)
  }

  out <- id_mat
  nz <- which(!is.na(id_mat) & id_mat > 0)
  out[nz] <- vapply(id_mat[nz], collapse_one, integer(1))
  out
}

# ----------------------------
# Matrix utilities
# ----------------------------

extract_coronal_slice <- function(vol,
                                  section,
                                  reverse_ap = FALSE,
                                  reverse_dv = FALSE,
                                  reverse_ml = FALSE) {
  d <- dim(vol)
  if (length(d) != 3L) stop("Expected a 3D volume")
  if (section < 1L || section > d[1]) stop("Section out of bounds: ", section, " not in 1..", d[1])

  sec <- if (reverse_ap) d[1] - section + 1L else section
  mat <- vol[sec, , , drop = TRUE]  # DV x ML
  if (reverse_dv) mat <- mat[nrow(mat):1L, , drop = FALSE]
  if (reverse_ml) mat <- mat[, ncol(mat):1L, drop = FALSE]
  mat
}

crop_to_mask <- function(mask, pad = 1L) {
  idx <- which(mask, arr.ind = TRUE)
  if (nrow(idx) == 0) return(list(r = c(1L, nrow(mask)), c = c(1L, ncol(mask))))
  rmin <- max(1L, min(idx[, 1]) - pad)
  rmax <- min(nrow(mask), max(idx[, 1]) + pad)
  cmin <- max(1L, min(idx[, 2]) - pad)
  cmax <- min(ncol(mask), max(idx[, 2]) + pad)
  list(r = c(rmin, rmax), c = c(cmin, cmax))
}

upsample_nearest <- function(mat, factor = 1L) {
  factor <- as.integer(factor)
  if (factor <= 1L) return(mat)
  rr <- rep(seq_len(nrow(mat)), each = factor)
  cc <- rep(seq_len(ncol(mat)), each = factor)
  mat[rr, cc, drop = FALSE]
}

upsample_bilinear <- function(mat, factor = 1L) {
  factor <- as.integer(factor)
  if (factor <= 1L) return(mat)
  nr <- nrow(mat)
  nc <- ncol(mat)
  x_old <- seq_len(nc)
  y_old <- seq_len(nr)
  x_new <- seq(1, nc, length.out = nc * factor)
  y_new <- seq(1, nr, length.out = nr * factor)

  m2 <- mat
  m2[!is.finite(m2)] <- NA_real_

  # fill rows with linear interpolation where possible
  tmp <- matrix(NA_real_, nrow = nr, ncol = length(x_new))
  for (i in seq_len(nr)) {
    row <- m2[i, ]
    keep <- is.finite(row)
    if (sum(keep) == 0) next
    if (sum(keep) == 1) {
      tmp[i, ] <- row[keep][1]
    } else {
      tmp[i, ] <- approx(x_old[keep], row[keep], xout = x_new, rule = 2)$y
    }
  }

  out <- matrix(NA_real_, nrow = length(y_new), ncol = length(x_new))
  for (j in seq_len(ncol(tmp))) {
    col <- tmp[, j]
    keep <- is.finite(col)
    if (sum(keep) == 0) next
    if (sum(keep) == 1) {
      out[, j] <- col[keep][1]
    } else {
      out[, j] <- approx(y_old[keep], col[keep], xout = y_new, rule = 2)$y
    }
  }

  out
}

boundary_from_ids <- function(id_mat) {
  a <- id_mat
  a[is.na(a)] <- 0
  nr <- nrow(a)
  nc <- ncol(a)
  b <- matrix(FALSE, nrow = nr, ncol = nc)

  if (nr > 1) {
    diff_ud <- a[-1, , drop = FALSE] != a[-nr, , drop = FALSE]
    b[-1, ] <- b[-1, ] | diff_ud
    b[-nr, ] <- b[-nr, ] | diff_ud
  }
  if (nc > 1) {
    diff_lr <- a[, -1, drop = FALSE] != a[, -nc, drop = FALSE]
    b[, -1] <- b[, -1] | diff_lr
    b[, -nc] <- b[, -nc] | diff_lr
  }
  b[a == 0] <- FALSE
  b
}

# ----------------------------
# Color utilities
# ----------------------------

to_expr_colors <- function(mat,
                           zlim = NULL,
                           palette = c("#f7f7ef", "#e5d48a", "#aa6b3f", "#000000"),
                           n = 256L) {
  pal <- grDevices::colorRampPalette(palette)(n)
  finite <- is.finite(mat)
  vals <- mat[finite]
  if (is.null(zlim)) {
    zlim <- range(vals, na.rm = TRUE, finite = TRUE)
  }
  if (!all(is.finite(zlim)) || zlim[1] == zlim[2]) {
    zlim <- c(0, 1)
  }
  idx <- round((mat - zlim[1]) / (zlim[2] - zlim[1]) * (n - 1L)) + 1L
  idx[idx < 1L] <- 1L
  idx[idx > n] <- n
  cols <- matrix(NA_character_, nrow = nrow(mat), ncol = ncol(mat))
  cols[finite] <- pal[idx[finite]]
  list(colors = cols, zlim = zlim, palette = pal)
}

# ----------------------------
# Expression volume
# ----------------------------

build_expression_volume <- function(genes,
                                    dataset_ids = NULL,
                                    plane = c("coronal", "sagittal"),
                                    metric = c("energy", "density", "intensity"),
                                    per_gene_transform = c("none", "zscore"),
                                    verbose = TRUE) {
  plane <- match.arg(plane)
  metric <- match.arg(metric)
  per_gene_transform <- match.arg(per_gene_transform)

  genes <- as.character(genes)

  if (is.null(dataset_ids)) {
    ds_list <- lapply(genes, function(g) {
      ds <- find_section_datasets(
        gene = g,
        plane = plane,
        orientation = "antisense",
        require_ish = TRUE,
        verbose = verbose
      )
      ds$requested_gene <- g
      ds
    })
    ds_tbl <- do.call(rbind, ds_list)

    if (nrow(ds_tbl) == 0) stop("No experiments found for requested genes")

    # one dataset per requested gene, prefer highest id if multiple remain
    picked <- do.call(rbind, lapply(split(ds_tbl, ds_tbl$requested_gene), function(x) {
      x <- x[order(x$id, decreasing = TRUE), , drop = FALSE]
      x[1, , drop = FALSE]
    }))
    picked <- picked[match(genes, picked$requested_gene), , drop = FALSE]
    dataset_ids <- picked$id

    if (verbose) {
      message("Selected Allen experiments:")
      for (i in seq_along(genes)) {
        message("  ", genes[i], " -> ", dataset_ids[i])
      }
    }
  } else {
    dataset_ids <- as.integer(dataset_ids)
    if (length(dataset_ids) != length(genes)) {
      stop("dataset_ids must have same length as genes")
    }
    if (verbose) {
      message("Using supplied dataset_ids:")
      for (i in seq_along(genes)) {
        message("  ", genes[i], " -> ", dataset_ids[i])
      }
    }
  }

  vols <- lapply(dataset_ids, function(id) download_expression_grid(id, metric = metric, verbose = verbose))
  dims <- dim(vols[[1]])
  same_dims <- all(vapply(vols, function(v) identical(dim(v), dims), logical(1)))
  if (!same_dims) stop("Grid dimensions differ across genes")

  if (per_gene_transform == "zscore") {
    vols <- lapply(vols, function(v) {
      m <- mean(v, na.rm = TRUE)
      s <- stats::sd(v, na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 1
      (v - m) / s
    })
  }

  volume <- Reduce(`+`, vols) / length(vols)

  list(
    volume = volume,
    dataset_ids = dataset_ids,
    genes = genes
  )
}

# ----------------------------
# Panel raster composition
# ----------------------------

compose_section_raster <- function(expr_slice,
                                   ann_slice,
                                   structure_df,
                                   atlas_depth = 5L,
                                   smooth_factor = 4L,
                                   expr_interp = c("nearest", "bilinear"),
                                   expr_palette = c("#f7f7ef", "#e5d48a", "#aa6b3f", "#000000"),
                                   zlim = NULL,
                                   background = "#E6E6E6",
                                   silhouette = "#F2F2F2",
                                   boundary_col = "#B8B0A5") {
  expr_interp <- match.arg(expr_interp)

  mask <- is.finite(ann_slice) & ann_slice > 0
  expr_slice[!mask] <- NA_real_
  ann_slice[!mask] <- NA

  bbox <- crop_to_mask(mask, pad = 1L)
  expr_slice <- expr_slice[bbox$r[1]:bbox$r[2], bbox$c[1]:bbox$c[2], drop = FALSE]
  ann_slice <- ann_slice[bbox$r[1]:bbox$r[2], bbox$c[1]:bbox$c[2], drop = FALSE]
  mask <- mask[bbox$r[1]:bbox$r[2], bbox$c[1]:bbox$c[2], drop = FALSE]

  ann_disp <- collapse_ids_to_depth(ann_slice, structure_df, atlas_depth = atlas_depth)

  nc <- ncol(ann_disp)
  mid <- floor(nc / 2)

  ann_left <- ann_disp
  if (mid < nc) ann_left[, (mid + 1L):nc] <- NA

  expr_right <- expr_slice
  if (mid >= 1L) expr_right[, seq_len(mid)] <- NA
  expr_right[!mask] <- NA_real_

  mask_up <- upsample_nearest(mask, factor = smooth_factor)
  ann_left_up <- upsample_nearest(ann_left, factor = smooth_factor)
  ann_full_up <- upsample_nearest(ann_disp, factor = smooth_factor)
  expr_up <- if (expr_interp == "nearest") {
    upsample_nearest(expr_right, factor = smooth_factor)
  } else {
    upsample_bilinear(expr_right, factor = smooth_factor)
  }

  # Reapply mask after interpolation to avoid rectangular artifacts
  expr_up[!mask_up] <- NA_real_

  color_map <- structure_df$color
  names(color_map) <- as.character(structure_df$id)
  atlas_cols <- matrix(NA_character_, nrow = nrow(ann_left_up), ncol = ncol(ann_left_up))
  keep_atlas <- !is.na(ann_left_up) & ann_left_up > 0
  atlas_cols[keep_atlas] <- color_map[as.character(ann_left_up[keep_atlas])]

  expr_cols <- to_expr_colors(expr_up, zlim = zlim, palette = expr_palette)
  bnd <- boundary_from_ids(ann_full_up)

  img <- matrix(background, nrow = nrow(ann_full_up), ncol = ncol(ann_full_up))
  img[mask_up] <- silhouette
  img[keep_atlas] <- atlas_cols[keep_atlas]
  keep_expr <- is.finite(expr_up)
  img[keep_expr] <- expr_cols$colors[keep_expr]
  img[bnd] <- boundary_col

  list(
    img = img,
    zlim = expr_cols$zlim
  )
}

draw_colorbar <- function(zlim,
                          palette = c("#f7f7ef", "#e5d48a", "#aa6b3f", "#000000"),
                          title = "energy",
                          n_ticks = 4L) {
  pal <- grDevices::colorRampPalette(palette)(256)

  # Use absolute margins in inches, not line-based mar, to avoid
  # 'figure margins too large' on narrow layout cells.
  op <- par(mai = c(0.18, 0.03, 0.28, 0.55), xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)

  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rasterImage(as.raster(matrix(rev(pal), ncol = 1)), 0.20, 0.06, 0.52, 0.94, interpolate = TRUE)
  rect(0.20, 0.06, 0.52, 0.94, border = "grey40", lwd = 0.5)
  text(0.36, 1.02, labels = title, xpd = NA, cex = 0.85, font = 2)

  ticks <- seq(zlim[1], zlim[2], length.out = n_ticks)
  y <- seq(0.06, 0.94, length.out = n_ticks)
  segments(0.56, y, 0.64, y, xpd = NA)
  text(0.70, y, labels = formatC(ticks, digits = 3, format = "fg"),
       adj = c(0, 0.5), cex = 0.78, xpd = NA)
}

# ----------------------------
# Quick sanity-check helper
# ----------------------------

test_allen_download <- function(dataset_id = 71924291L,
                                metric = c("energy", "density", "intensity"),
                                verbose = TRUE) {
  metric <- match.arg(metric)
  ann <- download_annotation_grid(verbose = verbose)
  expr <- download_expression_grid(dataset_id, metric = metric, verbose = verbose)
  list(
    annotation_dim = dim(ann),
    expression_dim = dim(expr),
    expression_range = range(expr, na.rm = TRUE, finite = TRUE),
    n_annotation_ids = length(unique(as.integer(ann[is.finite(ann) & ann > 0])))
  )
}

# Integrated linefix8 additions for allen_sections_grid_combined.R
# Adds / fixes:
#   1) gene signatures from list_genes, with missing genes skipped + reported
#   2) signature calculation = per-gene z-score, clip to [-clip, clip],
#      sum and divide by number of genes used
#   3) custom palettes for expression heatmap
#   4) multi-panel layout with user-specified n_rows / n_cols
#   5) robust handling when some genes are missing or some volumes contain NA

# ---- boundary drawing --------------------------------------------------------

draw_label_boundaries_segments <- function(label_mat,
                                           xleft = 0, xright = 1,
                                           ybottom = 0, ytop = 1,
                                           smooth_factor = 8,
                                           col = grDevices::adjustcolor("#7A7066", alpha.f = 0.55),
                                           lwd = 0.3) {
  if (is.null(label_mat) || !is.matrix(label_mat) || any(dim(label_mat) < 2)) {
    return(invisible(NULL))
  }

  sf <- max(1L, as.integer(round(smooth_factor)))
  M <- label_mat[rep(seq_len(nrow(label_mat)), each = sf),
                 rep(seq_len(ncol(label_mat)), each = sf),
                 drop = FALSE]
  M[is.na(M)] <- 0

  nr <- nrow(M)
  nc <- ncol(M)

  left_ids  <- M[, -nc, drop = FALSE]
  right_ids <- M[, -1,  drop = FALSE]
  vb <- (left_ids != right_ids) & !(left_ids == 0 & right_ids == 0)

  top_ids <- M[-nr, , drop = FALSE]
  bot_ids <- M[-1,  , drop = FALSE]
  hb <- (top_ids != bot_ids) & !(top_ids == 0 & bot_ids == 0)

  x_edges <- seq(xleft, xright, length.out = nc + 1L)
  y_edges <- seq(ytop, ybottom, length.out = nr + 1L)

  if (nc > 1L) {
    for (j in seq_len(nc - 1L)) {
      rows <- which(vb[, j])
      if (!length(rows)) next
      split_idx <- c(1L, which(diff(rows) > 1L) + 1L)
      starts <- rows[split_idx]
      ends <- c(rows[split_idx[-1L] - 1L], rows[length(rows)])
      x <- x_edges[j + 1L]
      segments(rep(x, length(starts)), y_edges[starts],
               rep(x, length(ends)),   y_edges[ends + 1L],
               col = col, lwd = lwd, lend = 1)
    }
  }

  if (nr > 1L) {
    for (i in seq_len(nr - 1L)) {
      cols <- which(hb[i, ])
      if (!length(cols)) next
      split_idx <- c(1L, which(diff(cols) > 1L) + 1L)
      starts <- cols[split_idx]
      ends <- c(cols[split_idx[-1L] - 1L], cols[length(cols)])
      y <- y_edges[i + 1L]
      segments(x_edges[starts],     rep(y, length(starts)),
               x_edges[ends + 1L],  rep(y, length(ends)),
               col = col, lwd = lwd, lend = 1)
    }
  }

  invisible(NULL)
}

# ---- helpers -----------------------------------------------------------------

unique_keep_order <- function(x) x[!duplicated(x)]

resolve_expr_palette <- function(expr_palette = c("#f7f7ef", "#e5d48a", "#aa6b3f", "#000000"),
                                 n = 7L) {
  default_pal <- c("#f7f7ef", "#e5d48a", "#aa6b3f", "#000000")

  if (is.null(expr_palette)) return(default_pal)

  if (is.function(expr_palette)) {
    out <- expr_palette(n)
    return(as.character(out))
  }

  expr_palette <- as.character(expr_palette)

  if (length(expr_palette) == 1L) {
    nm <- expr_palette[[1]]
    if (!nzchar(nm) || tolower(nm) == "allen") return(default_pal)

    hcl_ok <- tryCatch(nm %in% grDevices::hcl.pals(), error = function(e) FALSE)
    if (hcl_ok) return(grDevices::hcl.colors(n, palette = nm))

    col_ok <- tryCatch({ grDevices::col2rgb(nm); TRUE }, error = function(e) FALSE)
    if (col_ok) return(c("white", nm))

    warning("Unknown expr_palette preset '", nm, "'. Falling back to Allen palette.")
    return(default_pal)
  }

  expr_palette
}

resolve_section_layout <- function(n_panels,
                                   n_rows = NULL,
                                   n_cols = NULL,
                                   legend_panel_index = NULL,
                                   fill_by_row = TRUE) {
  if (is.null(n_rows) && is.null(n_cols)) {
    n_rows <- 1L
    n_cols <- n_panels
  } else if (is.null(n_rows)) {
    n_cols <- as.integer(n_cols)
    n_rows <- ceiling(n_panels / n_cols)
  } else if (is.null(n_cols)) {
    n_rows <- as.integer(n_rows)
    n_cols <- ceiling(n_panels / n_rows)
  } else {
    n_rows <- as.integer(n_rows)
    n_cols <- as.integer(n_cols)
  }

  if (!is.finite(n_rows) || !is.finite(n_cols) || n_rows < 1L || n_cols < 1L) {
    stop("n_rows and n_cols must be positive integers")
  }
  if (n_rows * n_cols < n_panels) {
    stop("n_rows * n_cols is smaller than the number of sections")
  }

  vals <- c(seq_len(n_panels), rep(0L, n_rows * n_cols - n_panels))
  lay_panels <- matrix(vals, nrow = n_rows, ncol = n_cols, byrow = fill_by_row)

  if (is.null(legend_panel_index)) legend_panel_index <- n_panels + 1L
  lay <- cbind(lay_panels, rep(legend_panel_index, n_rows))
  list(layout = lay, n_rows = n_rows, n_cols = n_cols, legend_index = legend_panel_index)
}

resolve_signature_inputs <- function(genes,
                                     dataset_ids = NULL,
                                     plane = "coronal",
                                     verbose = TRUE) {
  genes <- as.character(genes)
  genes <- trimws(genes)
  genes <- genes[nzchar(genes) & !is.na(genes)]
  genes <- unique_keep_order(genes)

  if (!length(genes)) stop("No genes supplied")

  used_genes <- character(0)
  used_ids <- integer(0)
  skipped_genes <- character(0)
  skipped_reason <- character(0)

  if (is.null(dataset_ids)) {
    for (g in genes) {
      ds <- suppressWarnings(tryCatch(
        find_section_datasets(
          gene = g,
          plane = plane,
          orientation = "antisense",
          require_ish = TRUE,
          verbose = FALSE
        ),
        error = function(e) NULL
      ))

      if (is.null(ds) || nrow(ds) == 0) {
        skipped_genes <- c(skipped_genes, g)
        skipped_reason <- c(skipped_reason, "not found in Allen atlas")
        if (verbose) message("Skipping gene '", g, "': not found in Allen atlas (plane = ", plane, ")")
        next
      }

      ds <- ds[order(ds$id, decreasing = TRUE), , drop = FALSE]
      used_genes <- c(used_genes, g)
      used_ids <- c(used_ids, as.integer(ds$id[1]))
    }
  } else {
    dataset_ids <- as.integer(dataset_ids)
    if (length(dataset_ids) != length(genes)) {
      stop("dataset_ids must have the same length as genes/list_genes")
    }
    for (i in seq_along(genes)) {
      g <- genes[i]
      id <- dataset_ids[i]
      if (!is.finite(id) || is.na(id) || id <= 0) {
        skipped_genes <- c(skipped_genes, g)
        skipped_reason <- c(skipped_reason, "missing dataset_id")
        if (verbose) message("Skipping gene '", g, "': missing/invalid dataset_id")
        next
      }
      used_genes <- c(used_genes, g)
      used_ids <- c(used_ids, id)
    }
  }

  if (!length(used_genes)) {
    stop("None of the supplied genes could be resolved in Allen atlas")
  }

  list(
    genes = used_genes,
    dataset_ids = used_ids,
    skipped_genes = skipped_genes,
    skipped_reason = skipped_reason
  )
}

build_expression_volume_signature_skipmissing <- function(genes,
                                                          dataset_ids = NULL,
                                                          plane = c("coronal", "sagittal"),
                                                          metric = c("energy", "density", "intensity"),
                                                          per_gene_transform = c("none", "zscore"),
                                                          signature_mode = c("auto", "zscore_clip_mean", "raw_mean"),
                                                          signature_clip = 3,
                                                          verbose = TRUE) {
  plane <- match.arg(plane)
  metric <- match.arg(metric)
  per_gene_transform <- match.arg(per_gene_transform)
  signature_mode <- match.arg(signature_mode)

  sig <- resolve_signature_inputs(
    genes = genes,
    dataset_ids = dataset_ids,
    plane = plane,
    verbose = verbose
  )

  vols <- list()
  keep_genes <- character(0)
  keep_ids <- integer(0)
  ref_dims <- NULL
  skipped_genes <- sig$skipped_genes
  skipped_reason <- sig$skipped_reason

  for (i in seq_along(sig$dataset_ids)) {
    g <- sig$genes[i]
    id <- sig$dataset_ids[i]

    v <- tryCatch(
      download_expression_grid(id, metric = metric, verbose = verbose),
      error = function(e) e
    )

    if (inherits(v, "error")) {
      skipped_genes <- c(skipped_genes, g)
      skipped_reason <- c(skipped_reason, paste0("grid download failed: ", conditionMessage(v)))
      if (verbose) message("Skipping gene '", g, "' (dataset ", id, "): ", conditionMessage(v))
      next
    }

    if (is.null(ref_dims)) {
      ref_dims <- dim(v)
    } else if (!identical(dim(v), ref_dims)) {
      skipped_genes <- c(skipped_genes, g)
      skipped_reason <- c(skipped_reason, "grid dimensions differ from other genes")
      if (verbose) message("Skipping gene '", g, "' (dataset ", id, "): grid dimensions differ")
      next
    }

    vols[[length(vols) + 1L]] <- v
    keep_genes <- c(keep_genes, g)
    keep_ids <- c(keep_ids, id)
  }

  if (!length(vols)) {
    stop("No valid Allen expression volumes remained after filtering/skipping")
  }

  use_signature <- switch(
    signature_mode,
    auto = length(vols) > 1L,
    zscore_clip_mean = TRUE,
    raw_mean = FALSE
  )

  if (use_signature) {
    # Requested signature logic:
    # 1) per-gene z-score
    # 2) values < -clip are set to -clip; values > clip are set to clip
    # 3) clipped z-scores are summed
    # 4) divided by the number of genes used
    # IMPORTANT FIX:
    #   if sd == 0 (all values identical, e.g. all zero), use z = 0 rather than NA
    #   so the region remains plotted at the lowest heatmap value.
    zvols <- lapply(vols, function(v) {
      m <- mean(v, na.rm = TRUE)
      s <- stats::sd(v, na.rm = TRUE)
      z <- array(NA_real_, dim = dim(v))
      ok <- is.finite(v)

      if (!is.finite(m) || !is.finite(s) || s == 0) {
        z[ok] <- 0
      } else {
        z[ok] <- (v[ok] - m) / s
        z[z < -signature_clip] <- -signature_clip
        z[z >  signature_clip] <-  signature_clip
      }
      z
    })

    denom <- length(zvols)
    sum_vol <- array(0, dim = ref_dims)
    finite_any <- array(FALSE, dim = ref_dims)
    for (z in zvols) {
      zz <- z
      finite_any <- finite_any | is.finite(zz)
      zz[!is.finite(zz)] <- 0
      sum_vol <- sum_vol + zz
    }

    volume <- sum_vol / denom
    # Only keep NA where absolutely no finite values were ever available.
    # Constant-zero genes remain 0 and are therefore plotted.
    volume[!finite_any] <- NA_real_
    mode_used <- "zscore_clip_mean"
  } else {
    if (per_gene_transform == "zscore") {
      vols <- lapply(vols, function(v) {
        m <- mean(v, na.rm = TRUE)
        s <- stats::sd(v, na.rm = TRUE)
        if (!is.finite(s) || s == 0) s <- 1
        z <- (v - m) / s
        z[!is.finite(z)] <- NA_real_
        z
      })
    }

    # average while ignoring NA to avoid propagating missing values
    denom_arr <- array(0, dim = ref_dims)
    sum_arr <- array(0, dim = ref_dims)
    for (v in vols) {
      ok <- is.finite(v)
      vv <- v
      vv[!ok] <- 0
      sum_arr <- sum_arr + vv
      denom_arr <- denom_arr + ok
    }
    volume <- sum_arr / pmax(denom_arr, 1)
    volume[denom_arr == 0] <- NA_real_
    mode_used <- if (length(vols) > 1L) "raw_mean" else per_gene_transform
  }

  if (verbose && length(skipped_genes)) {
    message("Skipped genes: ", paste(skipped_genes, collapse = ", "))
  }

  list(
    volume = volume,
    dataset_ids = keep_ids,
    genes = keep_genes,
    skipped_genes = skipped_genes,
    skipped_reason = skipped_reason,
    signature_mode = mode_used,
    signature_clip = signature_clip
  )
}

compose_section_raster_linefix <- function(expr_slice,
                                           ann_slice,
                                           structure_df,
                                           atlas_depth = 5L,
                                           smooth_factor = 4L,
                                           expr_interp = c("nearest", "bilinear"),
                                           expr_palette = c("#f7f7ef", "#e5d48a", "#aa6b3f", "#000000"),
                                           zlim = NULL,
                                           background = "#FFFFFF",
                                           silhouette = "#F2F2F2") {
  expr_interp <- match.arg(expr_interp)

  mask <- is.finite(ann_slice) & ann_slice > 0
  expr_slice[!mask] <- NA_real_
  ann_slice[!mask] <- NA

  bbox <- crop_to_mask(mask, pad = 1L)
  expr_slice <- expr_slice[bbox$r[1]:bbox$r[2], bbox$c[1]:bbox$c[2], drop = FALSE]
  ann_slice  <- ann_slice[bbox$r[1]:bbox$r[2], bbox$c[1]:bbox$c[2], drop = FALSE]
  mask       <- mask[bbox$r[1]:bbox$r[2], bbox$c[1]:bbox$c[2], drop = FALSE]

  ann_disp <- collapse_ids_to_depth(ann_slice, structure_df, atlas_depth = atlas_depth)

  nc <- ncol(ann_disp)
  mid <- floor(nc / 2)

  ann_left <- ann_disp
  if (mid < nc) ann_left[, (mid + 1L):nc] <- NA

  expr_right <- expr_slice
  if (mid >= 1L) expr_right[, seq_len(mid)] <- NA
  expr_right[!mask] <- NA_real_

  mask_up <- upsample_nearest(mask, factor = smooth_factor)
  ann_left_up <- upsample_nearest(ann_left, factor = smooth_factor)
  expr_up <- if (expr_interp == "nearest") {
    upsample_nearest(expr_right, factor = smooth_factor)
  } else {
    upsample_bilinear(expr_right, factor = smooth_factor)
  }
  expr_up[!mask_up] <- NA_real_

  color_map <- structure_df$color
  names(color_map) <- as.character(structure_df$id)
  atlas_cols <- matrix(NA_character_, nrow = nrow(ann_left_up), ncol = ncol(ann_left_up))
  keep_atlas <- !is.na(ann_left_up) & ann_left_up > 0
  atlas_cols[keep_atlas] <- color_map[as.character(ann_left_up[keep_atlas])]

  expr_cols <- to_expr_colors(expr_up, zlim = zlim, palette = expr_palette)

  img <- matrix(background, nrow = nrow(ann_left_up), ncol = ncol(ann_left_up))
  img[mask_up] <- silhouette
  img[keep_atlas] <- atlas_cols[keep_atlas]
  keep_expr <- is.finite(expr_up)
  img[keep_expr] <- expr_cols$colors[keep_expr]

  list(
    img = img,
    zlim = expr_cols$zlim,
    ann_disp = ann_disp
  )
}

# ---- main plotting function --------------------------------------------------

plot_allen_sections_grid <- function(genes,
                                     dataset_ids = NULL,
                                     sections = 22:24,
                                     plane = c("coronal", "sagittal"),
                                     metric = c("energy", "density", "intensity"),
                                     per_gene_transform = c("none", "zscore"),
                                     signature_mode = c("auto", "zscore_clip_mean", "raw_mean"),
                                     signature_clip = 3,
                                     atlas_depth = 5L,
                                     smooth_factor = 4L,
                                     expr_interp = c("nearest", "bilinear"),
                                     reverse_ap = FALSE,
                                     reverse_dv = FALSE,
                                     reverse_ml = FALSE,
                                     expr_palette = c("#f7f7ef", "#e5d48a", "#aa6b3f", "#000000"),
                                     out_file = NULL,
                                     width = 12,
                                     height = 3.8,
                                     res = 220,
                                     legend_panel_width = 0.28,
                                     verbose = TRUE,
                                     boundary_lwd = 0.3,
                                     boundary_col = grDevices::adjustcolor("#7A7066", alpha.f = 0.55),
                                     main_title = NULL,
                                     show_plot = TRUE,
                                     list_genes = NULL,
                                     skip_missing_genes = TRUE,
                                     n_rows = NULL,
                                     n_cols = NULL,
                                     fill_by_row = TRUE) {
  plane <- match.arg(plane)
  metric <- match.arg(metric)
  per_gene_transform <- match.arg(per_gene_transform)
  expr_interp <- match.arg(expr_interp)
  signature_mode <- match.arg(signature_mode)

  if (plane != "coronal") {
    stop("This patch currently implements coronal plotting only")
  }

  if (!is.null(list_genes)) {
    genes <- unlist(list_genes, use.names = FALSE)
  }
  genes <- as.character(genes)
  genes <- genes[nzchar(trimws(genes))]
  if (!length(genes)) stop("No genes supplied")

  expr_palette <- resolve_expr_palette(expr_palette)

  ann_vol <- download_annotation_grid(verbose = verbose)
  structure_df <- download_structure_graph(graph_id = 1L, verbose = verbose)

  expr <- build_expression_volume_signature_skipmissing(
    genes = genes,
    dataset_ids = dataset_ids,
    plane = plane,
    metric = metric,
    per_gene_transform = per_gene_transform,
    signature_mode = signature_mode,
    signature_clip = signature_clip,
    verbose = verbose
  )
  expr_vol <- expr$volume

  if (!identical(dim(ann_vol), dim(expr_vol))) {
    stop("annotation and expression volumes have different dimensions: ",
         paste(dim(ann_vol), collapse = "x"), " vs ", paste(dim(expr_vol), collapse = "x"))
  }

  ann_slices <- lapply(sections, function(sec) {
    extract_coronal_slice(ann_vol, sec,
                          reverse_ap = reverse_ap,
                          reverse_dv = reverse_dv,
                          reverse_ml = reverse_ml)
  })
  expr_slices <- lapply(sections, function(sec) {
    extract_coronal_slice(expr_vol, sec,
                          reverse_ap = reverse_ap,
                          reverse_dv = reverse_dv,
                          reverse_ml = reverse_ml)
  })

  vals <- c()
  for (i in seq_along(sections)) {
    mask <- is.finite(ann_slices[[i]]) & ann_slices[[i]] > 0
    x <- expr_slices[[i]]
    x[!mask] <- NA_real_
    nc <- ncol(x)
    mid <- floor(nc / 2)
    if (mid >= 1L) x[, seq_len(mid)] <- NA_real_
    vals <- c(vals, x[is.finite(x)])
  }
  if (!length(vals)) {
    # fallback to all finite values from full expression volume
    vals <- expr_vol[is.finite(expr_vol)]
  }
  zlim <- range(vals, na.rm = TRUE, finite = TRUE)
  if (!length(vals) || !all(is.finite(zlim)) || zlim[1] == zlim[2]) zlim <- c(0, 1)

  panels <- lapply(seq_along(sections), function(i) {
    compose_section_raster_linefix(
      expr_slice = expr_slices[[i]],
      ann_slice = ann_slices[[i]],
      structure_df = structure_df,
      atlas_depth = atlas_depth,
      smooth_factor = smooth_factor,
      expr_interp = expr_interp,
      expr_palette = expr_palette,
      zlim = zlim
    )
  })

  panel_layout <- resolve_section_layout(
    n_panels = length(sections),
    n_rows = n_rows,
    n_cols = n_cols,
    fill_by_row = fill_by_row
  )

  draw_once <- function() {
    graphics::layout(panel_layout$layout,
                     widths = c(rep(1, panel_layout$n_cols), legend_panel_width),
                     heights = rep(1, panel_layout$n_rows))
    op <- par(bg = "#E6E6E6")
    on.exit(par(op), add = TRUE)

    if (!is.null(main_title) && nzchar(main_title)) {
      par(oma = c(0, 0, 2.2, 0))
    }

    for (i in seq_along(sections)) {
      par(mai = c(0.02, 0.02, 0.28, 0.02), xaxs = "i", yaxs = "i")
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
      rasterImage(as.raster(panels[[i]]$img), 0, 0, 1, 1, interpolate = FALSE)
      draw_label_boundaries_segments(panels[[i]]$ann_disp,
                                     xleft = 0, xright = 1,
                                     ybottom = 0, ytop = 1,
                                     smooth_factor = smooth_factor,
                                     col = boundary_col,
                                     lwd = boundary_lwd)
      title(main = as.character(sections[i]),
            col.main = "#E15759", font.main = 2, cex.main = 1.35, line = 0.05)
    }

    draw_colorbar(zlim = zlim, palette = expr_palette, title = metric)

    if (!is.null(main_title) && nzchar(main_title)) {
      mtext(main_title, outer = TRUE, side = 3, line = 0.4, font = 2, cex = 1.25)
    }
  }

  if (!is.null(out_file)) {
    grDevices::png(filename = out_file,
                   width = width, height = height,
                   units = "in", res = res, bg = "white")
    draw_once()
    grDevices::dev.off()
    if (show_plot) draw_once()
  } else if (show_plot) {
    draw_once()
  }

  invisible(list(
    sections = sections,
    genes_requested = unique_keep_order(as.character(if (!is.null(list_genes)) unlist(list_genes, use.names = FALSE) else genes)),
    genes_used = expr$genes,
    genes_skipped = expr$skipped_genes,
    skipped_reason = expr$skipped_reason,
    dataset_ids = expr$dataset_ids,
    zlim = zlim,
    annotation_dim = dim(ann_vol),
    expression_dim = dim(expr_vol),
    signature_mode = expr$signature_mode,
    signature_clip = expr$signature_clip,
    layout = list(n_rows = panel_layout$n_rows, n_cols = panel_layout$n_cols)
  ))
}

message("allen_sections_grid.R loaded.")

