cd your-repo-name
mkdir -p "week 1"
cat > "week 1/KEYWORDS.md" <<'EOF'
# Key Topics — Week 1

**موضوع کلی:** مقالات منتخب در حوزه‌ی *شبکه‌های مخابراتی و شبکه‌های کامپیوتری* با کاربردهای **یادگیری ماشین (ML)** و **هوش مصنوعی (AI)** — مقالات 2023 به بعد، مجلات با Impact Factor > 2.

---

## ۵ کلمه‌کلیدی مهم

1. **Machine Learning (یادگیری ماشین)**
   - کاربرد: پیش‌بینی رفتار شبکه، تشخیص ناهنجاری، تصمیم‌گیری خودکار.
   - نمونه کاربرد در مقالات: supervised/unsupervised learning برای تحلیل ترافیک و تشخیص حملات.

2. **Wireless Communications (ارتباطات بی‌سیم)**
   - کاربرد: بهینه‌سازی لایه‌های فیزیکی و لینک در 5G/6G، مدیریت طیف و QoS.
   - نمونه کاربرد: تخصیص منابع در شبکه‌های سلولی و IoT.

3. **Resource Allocation (تخصیص منابع)**
   - کاربرد: بهینه‌سازی پهنای‌باند، توان، و اسلات‌های زمانی با استفاده از RL و الگوریتم‌های مبتنی بر ML.
   - نمونه کاربرد: یادگیری تقویتی برای scheduling و کنترل قدرت.

4. **Deep Learning (یادگیری عمیق)**
   - کاربرد: برآورد کانال، تشخیص الگوهای پیچیده، طراحی گیرنده‌های end-to-end.
   - نمونه کاربرد: شبکه‌های عصبی عمیق برای channel estimation و signal detection.

5. **Edge Intelligence / Edge Computing (هوش لبه‌ای / محاسبات لبه‌ای)**
   - کاربرد: انجام inference و آموزش جزئی در لبه برای کاهش تأخیر و مصرف پهنای‌باند.
   - نمونه کاربرد: federated learning و distributed inference در نودهای لبه.

---

## تگ‌ها (برای استفاده در GitHub / Repo)
`machine-learning`, `deep-learning`, `wireless-communications`, `resource-allocation`, `edge-computing`, `networking`, `telecommunications`

---

## یادداشت
- این فایل می‌تواند به‌عنوان مرجع سریع (summary) برای پوشه‌ی `week 1` استفاده شود.
- در صورت نیاز می‌تونم همین فایل رو با فرمت فارسی کامل‌تر یا انگلیسی با فرمت IEEE-bibliography به‌روز کنم.
EOF

git add "week 1/KEYWORDS.md"
git commit -m "Add KEYWORDS.md for Week 1 (ML/AI in telecom & networking)"
git push
