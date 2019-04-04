using Documenter, wellwisejl

makedocs(
    modules = [wellwisejl],
    format = :html,
    checkdocs = :exports,
    sitename = "wellwisejl",
    pages = Any["index.md"],
)

deploydocs(
    repo = "github.com/alexpkeil1/wellwisejl.git",
)