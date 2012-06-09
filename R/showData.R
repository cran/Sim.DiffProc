# Sun May 13 13:30:43 2012
# Arsalane Chouaib GUIDOUM <starsalane@gmail.com>

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

showData <-
function (dataframe, colname.bgcolor = "grey50", rowname.bgcolor = "grey50", 
    body.bgcolor = "white", colname.textcolor = "white", rowname.textcolor = "white", 
    body.textcolor = "black", font = "Courier 12", maxheight = 30, 
    maxwidth = 80, title = NULL, rowname.bar = "left", colname.bar = "top", 
    rownumbers = FALSE, placement = "-20-40", suppress.X11.warnings = TRUE) 
{
    object.name <- deparse(substitute(dataframe))
    if (!is.data.frame(dataframe)) {
        temp <- try(dataframe <- as.data.frame(dataframe), silent = FALSE)
        if (inherits(temp, "try-error")) {
            stop(paste(object.name, "cannot be coerced to a data frame"))
        }
        object.name <- paste("as.data.frame(", object.name, ")", 
            sep = "")
    }
    if (is.numeric(rownumbers) && length(rownumbers) != nrow(dataframe)) {
        stop("rownumbers argument must be TRUE, FALSE or have length nrow(dataframe)")
    }
    oldwidth <- unlist(options("width"))
    options(width = 10000)
    conn <- file()
    sink(conn)
    print(dataframe)
    sink()
    zz <- scan(conn, sep = "\n", what = character(0), quiet = TRUE)
    close(conn)
    if (length(zz) > 1 + nrow(dataframe)) 
        stop("data frame too wide")
    options(width = oldwidth)
    if (suppress.X11.warnings) {
        messages.connection <- textConnection(".messages", open = "w", 
            local = TRUE)
        sink(messages.connection, type = "message")
        on.exit({
            sink(type = "message")
            close(messages.connection)
        })
    }
	tclRequire("BWidget")
    tclRequire("Tktable")
    base <- tktoplevel()
    tkwm.iconbitmap(base, default = file.path(R.home("bin"), "Rgui.exe"))
    tkwm.geometry(base, placement)
    tkwm.title(base, {
        if (is.null(title)) 
            object.name
        else title
    })
    nrows <- length(zz) - 1
    if (is.numeric(rownumbers)) 
        rowname.text <- paste(rownumbers, row.names(dataframe))
    else if (rownumbers) 
        rowname.text <- paste(1:nrows, row.names(dataframe))
    else rowname.text <- row.names(dataframe)
    namewidth = max(nchar(rowname.text))
    yy <- substring(zz, 2 + max(nchar(row.names(dataframe))))
    datawidth <- max(nchar(yy))
    winwidth <- min(1 + datawidth, maxwidth)
    hdr <- tktext(base, bg = colname.bgcolor, fg = colname.textcolor, 
        font = font, height = 1, width = winwidth, takefocus = TRUE)
    ftr <- tktext(base, bg = colname.bgcolor, fg = colname.textcolor, 
        font = font, height = 1, width = winwidth, takefocus = TRUE)
    textheight <- min(maxheight, nrows)
    txt <- tktext(base, bg = body.bgcolor, fg = body.textcolor, 
        font = font, height = textheight, width = winwidth, setgrid = 1, 
        takefocus = TRUE)
    lnames <- tktext(base, bg = rowname.bgcolor, fg = rowname.textcolor, 
        font = font, height = textheight, width = namewidth, 
        takefocus = TRUE)
    rnames <- tktext(base, bg = rowname.bgcolor, fg = rowname.textcolor, 
        font = font, height = textheight, width = namewidth, 
        takefocus = TRUE)
    xscroll <- tkscrollbar(base, orient = "horizontal", repeatinterval = 1, 
        command = function(...) {
            tkxview(txt, ...)
            tkxview(hdr, ...)
            tkxview(ftr, ...)
        })
    string.to.vector <- function(string.of.indices) {
        string.of.indices <- tclvalue(string.of.indices)
        as.numeric(strsplit(string.of.indices, split = " ")[[1]])
    }
    tkconfigure(txt, xscrollcommand = function(...) {
        tkset(xscroll, ...)
        xy <- string.to.vector(tkget(xscroll))
        tkxview.moveto(hdr, xy[1])
        tkxview.moveto(ftr, xy[1])
    })
    tkconfigure(hdr, xscrollcommand = function(...) {
        tkset(xscroll, ...)
        xy <- string.to.vector(tkget(xscroll))
        tkxview.moveto(txt, xy[1])
        tkxview.moveto(ftr, xy[1])
    })
    tkconfigure(ftr, xscrollcommand = function(...) {
        tkset(xscroll, ...)
        xy <- string.to.vector(tkget(xscroll))
        tkxview.moveto(hdr, xy[1])
        tkxview.moveto(txt, xy[1])
    })
    yscroll <- tkscrollbar(base, orient = "vertical", repeatinterval = 1, 
        command = function(...) {
            tkyview(txt, ...)
            tkyview(lnames, ...)
            tkyview(rnames, ...)
        })
    tkconfigure(txt, yscrollcommand = function(...) {
        tkset(yscroll, ...)
        xy <- string.to.vector(tkget(yscroll))
        tkyview.moveto(lnames, xy[1])
        tkyview.moveto(rnames, xy[1])
    })
    tkconfigure(lnames, yscrollcommand = function(...) {
        tkset(yscroll, ...)
        xy <- string.to.vector(tkget(yscroll))
        tkyview.moveto(txt, xy[1])
        tkyview.moveto(rnames, xy[1])
    })
    tkconfigure(rnames, yscrollcommand = function(...) {
        tkset(yscroll, ...)
        xy <- string.to.vector(tkget(yscroll))
        tkyview.moveto(txt, xy[1])
        tkyview.moveto(lnames, xy[1])
    })
    tkbind(txt, "<B2-Motion>", function(x, y) {
        tkscan.dragto(txt, x, y)
    })
    {
        copyText.hdr <- function() {
            tcl("event", "generate", .Tk.ID(hdr), "<<Copy>>")
        }
        tkbind(hdr, "<Button-1>", function() tkfocus(hdr))
        editPopupMenu.hdr <- tkmenu(hdr, tearoff = FALSE)
        tkadd(editPopupMenu.hdr, "command", label = "Copy <Ctrl-C>", 
            command = copyText.hdr)
        RightClick.hdr <- function(x, y) {
            rootx <- as.integer(tkwinfo("rootx", hdr))
            rooty <- as.integer(tkwinfo("rooty", hdr))
            xTxt <- as.integer(x) + rootx
            yTxt <- as.integer(y) + rooty
            tcl("tk_popup", editPopupMenu.hdr, xTxt, yTxt)
        }
        tkbind(hdr, "<Button-3>", RightClick.hdr)
        tkbind(hdr, "<Control-KeyPress-c>", copyText.hdr)
        copyText.ftr <- function() {
            tcl("event", "generate", .Tk.ID(ftr), "<<Copy>>")
        }
        tkbind(ftr, "<Button-1>", function() tkfocus(ftr))
        editPopupMenu.ftr <- tkmenu(ftr, tearoff = FALSE)
        tkadd(editPopupMenu.ftr, "command", label = "Copy <Ctrl-C>", 
            command = copyText.ftr)
        RightClick.ftr <- function(x, y) {
            rootx <- as.integer(tkwinfo("rootx", ftr))
            rooty <- as.integer(tkwinfo("rooty", ftr))
            xTxt <- as.integer(x) + rootx
            yTxt <- as.integer(y) + rooty
            tcl("tk_popup", editPopupMenu.ftr, xTxt, yTxt)
        }
        tkbind(ftr, "<Button-3>", RightClick.ftr)
        tkbind(ftr, "<Control-KeyPress-c>", copyText.ftr)
        copyText.txt <- function() {
            tcl("event", "generate", .Tk.ID(txt), "<<Copy>>")
        }
        tkbind(txt, "<Button-1>", function() tkfocus(txt))
        editPopupMenu.txt <- tkmenu(txt, tearoff = FALSE)
        tkadd(editPopupMenu.txt, "command", label = "Copy <Ctrl-C>", 
            command = copyText.txt)
        RightClick.txt <- function(x, y) {
            rootx <- as.integer(tkwinfo("rootx", txt))
            rooty <- as.integer(tkwinfo("rooty", txt))
            xTxt <- as.integer(x) + rootx
            yTxt <- as.integer(y) + rooty
            tcl("tk_popup", editPopupMenu.txt, xTxt, yTxt)
        }
        tkbind(txt, "<Button-3>", RightClick.txt)
        tkbind(txt, "<Control-KeyPress-c>", copyText.txt)
        copyText.lnames <- function() {
            tcl("event", "generate", .Tk.ID(lnames), "<<Copy>>")
        }
        tkbind(lnames, "<Button-1>", function() tkfocus(lnames))
        editPopupMenu.lnames <- tkmenu(lnames, tearoff = FALSE)
        tkadd(editPopupMenu.lnames, "command", label = "Copy <Ctrl-C>", 
            command = copyText.lnames)
        RightClick.lnames <- function(x, y) {
            rootx <- as.integer(tkwinfo("rootx", lnames))
            rooty <- as.integer(tkwinfo("rooty", lnames))
            xTxt <- as.integer(x) + rootx
            yTxt <- as.integer(y) + rooty
            tcl("tk_popup", editPopupMenu.lnames, xTxt, yTxt)
        }
        tkbind(lnames, "<Button-3>", RightClick.lnames)
        tkbind(lnames, "<Control-KeyPress-c>", copyText.lnames)
        copyText.rnames <- function() {
            tcl("event", "generate", .Tk.ID(rnames), "<<Copy>>")
        }
        tkbind(rnames, "<Button-1>", function() tkfocus(rnames))
        editPopupMenu.rnames <- tkmenu(rnames, tearoff = FALSE)
        tkadd(editPopupMenu.rnames, "command", label = "Copy <Ctrl-C>", 
            command = copyText.rnames)
        RightClick.rnames <- function(x, y) {
            rootx <- as.integer(tkwinfo("rootx", rnames))
            rooty <- as.integer(tkwinfo("rooty", rnames))
            xTxt <- as.integer(x) + rootx
            yTxt <- as.integer(y) + rooty
            tcl("tk_popup", editPopupMenu.rnames, xTxt, yTxt)
        }
        tkbind(rnames, "<Button-3>", RightClick.rnames)
        tkbind(rnames, "<Control-KeyPress-c>", copyText.rnames)
    }
    tktag.configure(hdr, "notwrapped", wrap = "none")
    tktag.configure(ftr, "notwrapped", wrap = "none")
    tktag.configure(txt, "notwrapped", wrap = "none")
    tktag.configure(lnames, "notwrapped", wrap = "none")
    tktag.configure(rnames, "notwrapped", wrap = "none")
    tkinsert(txt, "end", paste(paste(yy[-1], collapse = "\n"), 
        sep = ""), "notwrapped")
    tkgrid(txt, row = 1, column = 1, sticky = "nsew")
    if ("top" %in% colname.bar) {
        tkinsert(hdr, "end", paste(yy[1], sep = ""), "notwrapped")
        tkgrid(hdr, row = 0, column = 1, sticky = "ew")
    }
    if ("bottom" %in% colname.bar) {
        tkinsert(ftr, "end", paste(yy[1], sep = ""), "notwrapped")
        tkgrid(ftr, row = 2, column = 1, sticky = "ew")
    }
    if ("left" %in% rowname.bar) {
        tkinsert(lnames, "end", paste(rowname.text, collapse = "\n"), 
            "notwrapped")
        tkgrid(lnames, row = 1, column = 0, sticky = "ns")
    }
    if ("right" %in% rowname.bar) {
        tkinsert(rnames, "end", paste(rowname.text, collapse = "\n"), 
            "notwrapped")
        tkgrid(rnames, row = 1, column = 2, sticky = "ns")
    }
    tkconfigure(txt, state = "disabled")
    tkconfigure(lnames, state = "disabled")
    tkconfigure(rnames, state = "disabled")
    if (maxheight < nrows) {
        tkgrid(yscroll, row = 1, column = 3, sticky = "ns")
    }
    if (maxwidth < datawidth) {
        tkgrid(xscroll, row = 3, column = 1, sticky = "ew")
    }
    tkgrid.rowconfigure(base, 1, weight = 1)
    tkgrid.columnconfigure(base, 1, weight = 1)
    tkwm.maxsize(base, 1 + datawidth, nrows)
    tkwm.minsize(base, 1 + nchar(names(dataframe)[1]), 1)
    invisible(NULL)
}

