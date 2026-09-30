# Submitting `didc` to SSC

SSC (Statistical Software Components, Boston College) is the archive behind
`ssc install` and `findit`. Getting listed there is what makes a Stata command
discoverable to people who are not already reading the repository. The process
is an email, not a pull request, and there is no review of the statistics — but
the package must be self-contained and self-documenting.

Everything needed is already in this repository. This file is the checklist and
the text to send.

## What SSC requires

| requirement | status |
|---|---|
| a `.pkg` file describing the package | `didc.pkg` |
| one `.ado` file per command | `didc.ado`, `didc_test.ado`, `didc_bounds.ado` |
| internal helper `.ado` files | `_didc_prep.ado`, `_didc_lpoly.ado`, `_didc_ks.ado`, `_didc_display.ado` |
| a `.sthlp` help file for every documented command | `didc.sthlp`, `didc_test.sthlp`, `didc_bounds.sthlp` |
| no dependency that SSC will not install | depends on `rdrobust`, itself on SSC |
| a stated license | AGPL-3.0, in `LICENSE` |
| no illegal filenames (SSC wants lowercase, no spaces) | all files lowercase |

SSC hosts the files itself, so the `.ado` and `.sthlp` files must be sent as
attachments. The `data/` and `examples/` directories are **not** sent; SSC
archives code, not data, and the archive rejects large attachments.

Note that SSC puts files in a flat directory, which is why this repository uses
a flat layout for the `.ado` and `.sthlp` files already.

## How to submit

Send an email to `ssc@stata.com` with the attached files and roughly the
following text.

> Subject: `didc: Stata module for difference-in-discontinuities`
>
> Dear SSC,
>
> I would like to submit the package `didc` to the SSC archive.
>
> `didc` estimates, tests and partially identifies effects in a
> difference-in-discontinuities design: a regression-discontinuity setting where
> a confounding policy is assigned by the same cutoff as the treatment of
> interest, so that the continuity assumption underlying standard RD fails.
> Following Picchetti, Pinto & Shinoki (2026, arXiv:2405.18531), the estimator is
> a local polynomial regression of the change in the outcome over time at the
> cutoff; the package adds the paper's two validity tests and its
> partial-identification bounds.
>
> The package contains three user commands and four internal helpers, each with
> a help file. It depends only on `rdrobust`, which is on SSC.
>
> A repository with simulation data, an end-to-end example and a test suite that
> runs without R is at https://github.com/niuniuhaoyu/niuniuhaoyu-didc
>
> The package is licensed AGPL-3.0.
>
> Attached: didc.pkg, didc.ado, didc_test.ado, didc_bounds.ado,
> _didc_prep.ado, _didc_lpoly.ado, _didc_ks.ado, _didc_display.ado,
> didc.sthlp, didc_test.sthlp, didc_bounds.sthlp.
>
> Thank you,
> Haoyu Niu

## The `d` lines currently in `didc.pkg`

The archive's summary line is the first line of the `.pkg`. It is currently:

```
d 'didc': module for difference-in-discontinuities designs
```

The rest of the `.pkg` describes the three commands, the author, the
distribution date and the requirements. Keep the `d` lines under about 80
characters; the archive reformats long abstracts.

## Keywords to suggest

`didc`, `difference-in-discontinuities`, `regression discontinuity`,
`difference-in-differences`, `causal inference`, `policy evaluation`,
`partial identification`, `sensitivity analysis`.

## Before submitting

- [ ] the version in `didc.pkg`'s `Distribution-Date` matches the tag
- [ ] `CHANGELOG.md` has an entry for the version
- [ ] `net install didc, from(".../v0.1.0/")` has been tested from a clean
      machine (or a clean `plus` directory)
- [ ] every command has a `sthlp`, and `help didc` opens it
- [ ] `rdrobust` is declared in the `.pkg`'s `Requires:` line

## After SSC accepts it

The install instruction in the README becomes one line:

```stata
ssc install didc, replace
```

and `didc.pkg` should still be kept in the repository, because
`net install ... from(<github URL>)` remains the way to get the development
version.

SSC assigns a number of the form `S45xxxx`; it is worth putting it in the README
next to the citation, since that number is what people cite.
