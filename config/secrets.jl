try
  Genie.Util.isprecompiling() || Genie.Secrets.secret_token!("bfd7087c6f6267ae598e7105baf7d08c928473cf31705c45bde47726e0d73994")
catch ex
  @error "Failed to generate secrets file: $ex"
end
